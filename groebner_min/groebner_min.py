from sage.all import *

from .tree import *
from .timing import *

import itertools as it
import os
import sys
import subprocess
import time
import platform
import argparse


ESPRESSO_EXECUTABLE = os.environ.get("ESPRESSO_EXECUTABLE") or ("../Espresso.exe" if "microsoft" in platform.release() else "espresso")


def parse_polynomial(polynomial):
    negate = False
    monomials = []
    for monomial in polynomial.monomials():
        if monomial == 1:
            negate = True
            continue
        variables = [Value(var) for var in monomial.variables()]
        monomial = And(variables)
        monomials.append(monomial)
    polynomial = Xor(monomials)
    if negate:
        polynomial = Not(polynomial)
    return polynomial


def sample_polynomial(p):
    var = p.variables()
    n = len(var)
    values = []
    for i in range(1 << n):
        bitstring = [(i >> j) & 1 for j in range(n)]
        assignment = dict(zip(var, bitstring))
        value = p.subs(assignment)
        if value == 1:
            values.append(bitstring)
    return var, values

def read_bits_csv(path, bit_num):
    with open(f"{path}/bit_{bit_num}.csv") as f:
        data = f.readlines()

    with open(f"{path}/bits.csv") as f:
        bits_data = f.readlines()

    relevant_bits = [int(x.strip()) for x in bits_data[bit_num].split(",") if x.strip() != '']
    n = len(relevant_bits)

    values = []
    for line in data:
        x, y = line.split(",")
        x = int(x)
        y = int(y)

        bin_repr = []
        for bits in relevant_bits:
            bin_repr.append((x >> bits) & 1)

        values.append((tuple(bin_repr), y))

    return relevant_bits, values

def read_espresso_in(path):
    with open(f"{path}") as f:
        data = f.readlines()

    values = []
    seen_inputs = set()
    n = -1
    for line in data:
        line = line.strip()
        if line[0] == "." or line[0] == "#":
            if line.startswith(".i"):
                n = int(line.split(" ")[1])
            continue

        x, y = line.split(" ")
        y = int(y)

        assert n == -1 or n == len(x)
        n = len(x)

        bin_repr = []
        for bit in x:
            assert bit == "1" or bit == "0", "found don't care bit, this is currently not supported.s"
            bin_repr.append(int(bit))

        values.append((tuple(bin_repr), y))
        seen_inputs.add(tuple(bin_repr))

    for bits in it.product([0, 1], repeat=n):
        if bits not in seen_inputs:
            values.append((bits, 0))


    return n, values

def read_blackbox():
    n = 6
    blackbox = lambda x: x[1] & ((x[0] ^ x[5]) & (x[2]) | (x[3] & (x[4] ^ x[5])))
    values = []
    for bits in it.product([0, 1], repeat=n):
        values.append((bits, blackbox(bits)))
    return (n, values)


def compute_groebner_basis(bits, values, R, stats: RunStatistics):

    if len(values) == 0:
        return []

    stats.add_groebner_call()

    n = len(bits)
    if stats.verbose:
        print(f"Starting Espresso ({n} bits, {len(values)} truth-values).")
    espresso = subprocess.Popen(ESPRESSO_EXECUTABLE, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    write = lambda line: espresso.stdin.write(line.encode())
    write(f".i {n}\n")
    write(f".o 1\n")
    for bitstring in values:
        bin_repr = "".join(map(str, bitstring))
        write(f"{bin_repr} {1}\n")
    write(".e\n")


    time_start = time.time()
    out, err = espresso.communicate()
    time_end = time.time()

    stats.add_time_in_espresso((time_end - time_start))

    if err:
        print("Error:", err.decode())
        sys.exit(1)

    relations = []

    for line in out.decode().split("\n"):
        if not line or line.startswith("."):
            continue
        var, res = line.split(" ")
        term = 1
        for i in range(len(var)):
            if var[i] == "-":
                continue
            elif var[i] == '1':
                term *= bits[i]
            elif var[i] == '0':
                term *= bits[i] + 1
        relations.append(term)

    idempotency_relations = [bits[i]**2 + bits[i] for i in range(n)]
    relations.extend(idempotency_relations)


    I = R.ideal(relations)

    time_start = time.time()

    B = I.groebner_basis()
    B = [b for b in B if b not in idempotency_relations]
    time_end = time.time()
    stats.add_time_in_groebner((time_end - time_start))

    return B


def minimize_rec(x, bitstrings, R, stats, negate=False, depth=0, size=10**10):
    basis = compute_groebner_basis(x, bitstrings, R, stats)
    basis_size = sum(parse_polynomial(g).size() for g in basis)
    if basis_size > size:
        return
    result = []
    for element in basis:
        factors = []
        # This is needed if PolyBoRi is used, as factorization is not natively supported
        RR = PolynomialRing(GF(2), "x", stats.n)
        for fac, _ in factor(RR(element)): # TODO check _
            if stats.verbose:
                print("-" * depth, fac)
            if fac.degree() <= 1 or depth > 5: # better criterion // check for termination
                if negate:
                    factors.append(Not(parse_polynomial(fac)))
                else:
                    factors.append(parse_polynomial(fac))
            else:
                fac_size = parse_polynomial(fac).size()
                new_bits, new_values = sample_polynomial(fac + 1)
                min_fac = minimize_rec(new_bits, new_values, R, stats, not negate, depth + 1, 2 * fac_size)
                #seen.add(fac)

                if min_fac is None or fac_size < min_fac.size():
                    if negate:
                        factors.append(Not(parse_polynomial(fac)))
                    else:
                        factors.append(parse_polynomial(fac))
                else:
                    factors.append(min_fac)
        if negate:
            result.append(Or(factors))
        else:
            result.append(And(factors))
    if negate:
        return And(result)
    else:
        return Or(result)


# use argparse
# python groebner-min.py --input <path to espresso input> --output <output file> [--negate]

def main():
    time_start = time.time()

    argparser = argparse.ArgumentParser(description="Groebner Bases for Boolean Function Minimization")
    argparser.add_argument("input", type=str, help="Path to espresso input file (PLA format).")
    argparser.add_argument("--output", type=str, help="Store an extra output file with the generated formula and some statistics.")
    argparser.add_argument("--negate", action="store_true", help="Negate the output formula")
    argparser.add_argument("--stats", action="store_true", help="Print statistics about this run")
    argparser.add_argument("-v", action="store_true", help="Print more verbose information")
    argparser.add_argument("--pb", action="store_true", help="Use PolyBoRi instead of Sage for Groebner basis computation.")


    args = argparser.parse_args()
    # what happens on error?
    path = args.input
    negate = args.negate

    stats = RunStatistics()
    stats.verbose = args.v
    stats.print_stats = args.stats


    output = args.output


    n, data = read_espresso_in(path)
    stats.n = n

    check_value = 1 if not negate else 0
    values = [ value for value, y in data if y == check_value ]

    if args.pb:
        R = BooleanPolynomialRing(n, "x", order='lex')
    else:
        R = PolynomialRing(GF(2), "x", n)

    x = R.gens()

    tree = minimize_rec(x, values, R, stats, negate=negate)
    tree = modify(tree, simplify_all)
    tree = modify(tree, simplify_all)

    for leaf in tree.leafs():
        leaf.value = x.index(leaf.value)

    for leaf in tree.leafs():
        leaf.value = f"x[{leaf.value}]"

    time_tree_simplification_start = time.time()
    tree.traverse(n, stats)
    stats.add_time_tree_simpl(time.time() - time_tree_simplification_start)


    extract_xors = False

    if stats.verbose:
        print("Generating Code.")
    if not extract_xors:
        code = "lambda x: " + str(tree) + " & 1\n"
        tree_size = tree.size()
    #else:
    #    tmp_tree = tree.copy()
    #    tree_size = 0
    #    code = "def formula(x):\n"
    #    remap = {}
    #    for node in tmp_tree.find("^"):
    #        str_node = str(node)
    #        if str_node in remap:
    #            idx = remap[str_node]
    #        else:
    #            idx = remap[str_node] = len(remap)
    #            code += f"    y{idx} = " + str_node + "\n"
    #            tree_size += node.size()
    #        node.children = [LeafNode(f"y{idx}")]
    #    code += f"    return " + str(tmp_tree) + " & 1\n"
    #    tmp_tree = modify(tmp_tree, simplify_all)
    #    tree_size += tmp_tree.size()

    print(f"formula = {code}", end="")

    formula = eval(code)

    if stats.verbose:
        print("Testing validity...")
    for value, y in data:
        if formula(value) != y:
            print(f"Error: formula({value}) != {y}")
            sys.exit(1)

    #
    import os
    input_path = sys.argv[1]

    stats.TIME_TOTAL += time.time() - time_start

    if stats.print_stats:
        print("\nTiming: ")
        print(f"{stats.TIME_IN_ESPRESSO=:.3f}")
        print(f"{stats.TIME_IN_GROEBNER=:.3f}")
        print(f"{stats.TIME_IN_SETCOVER=:.3f}")
        print(f"{stats.TIME_TREE_SIMPLIFICATION=:.3f}")
        print(f"{stats.TIME_TOTAL=:.3f}")
        print("Stats: ")
        print(f"{stats.TOTAL_GROEBNER_CALLS=}")
        print()
        print(f"Number of input bits {n=}")
        print(f"Final formula size: {tree_size}")

    if output:
        with open(f"{dir}/{output}", "w") as f:
            f.write(code)
            f.write("Timing: \n")
            f.write(f"TIME_IN_ESPRESSO={stats.TIME_IN_ESPRESSO:.3f}\n")
            f.write(f"TIME_IN_GROEBNER={stats.TIME_IN_GROEBNER:.3f}\n")
            f.write(f"TIME_IN_SETCOVER={stats.TIME_IN_SETCOVER:.3f}\n")
            f.write(f"TIME_TREE_SIMPLIFICATION={stats.TIME_TREE_SIMPLIFICATION:.3f}\n")
            f.write(f"TIME_TOTAL={stats.TIME_TOTAL:.3f}\n")
            f.write("\nStats: \n")
            f.write(f"{stats.TOTAL_GROEBNER_CALLS=}\n")
            f.write(f"{n=}\n")
            f.write(f"{tree_size=}\n")
            f.write("Done.")

    if stats.verbose:
        print("Done.")

if __name__ == "__main__":
    main()