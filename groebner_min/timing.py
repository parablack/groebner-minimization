import time

class RunStatistics:
    def __init__(self):
        self.TIME_TOTAL = 0
        self.TIME_IN_ESPRESSO = 0
        self.TIME_IN_GROEBNER = 0
        self.TIME_IN_SETCOVER = 0
        self.TIME_TREE_SIMPLIFICATION = 0
        self.TOTAL_GROEBNER_CALLS = 0
        self.n = 0
        self.verbose = False
        self.print_stats = False

    def add_groebner_call(self):
        self.TOTAL_GROEBNER_CALLS += 1

    def add_time_in_espresso(self, time):
        self.TIME_IN_ESPRESSO += time

    def add_time_in_groebner(self, time):
        self.TIME_IN_GROEBNER += time

    def add_time_in_setcover(self, time):
        self.TIME_IN_SETCOVER += time

    def add_time_tree_simpl(self, time):
        self.TIME_TREE_SIMPLIFICATION += time
