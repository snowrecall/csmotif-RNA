from numpy import array

#===============================================================================
class BCTab:
    """Barcode -> chemical-shift table container.

    The table file is whitespace-delimited with a header row; the first
    column is a barcode (e.g. 'GC-AU-AU') and the following columns are
    numeric chemical shifts.  The shift type (NH, CH, ...) is inferred
    from the value ranges.
    """

    #---------------------------------------------------------------------------
    def __init__(self, tab=None):
        self.CH   = {}
        self.NH   = {}
        self.CHnu = {}
        self.NHnu = {}
        self.CHHz = {}
        self.NHHz = {}
        self.N    = {}
        self.H    = {}
        self.C    = {}
        self.num  = {}
        if tab:
            self.readTab(tab)

    #---------------------------------------------------------------------------
    def readTab(self, tab):
        # FIX (Py2 -> Py3): ``map`` returns an iterator in Python 3 and was
        # exhausted after the first comprehension, so Tab2 / self.num ended
        # up empty and self.H raised KeyError.  Materialise the list once.
        with open(tab) as fh:
            lines = [ln.strip() for ln in fh.readlines()[1:] if ln.strip()]

        Tab1 = {x.split()[0]: float(x.split()[1]) for x in lines}
        Tab2 = {x.split()[0]: float(x.split()[2]) for x in lines}
        try:
            self.num = {x.split()[0]: int(x.split()[3]) for x in lines}
        except IndexError:
            pass

        bds = sorted(Tab1.keys())
        v1 = array(list(Tab1.values()))
        v2 = array(list(Tab2.values()))

        # -------- 15N / 1H (imino) -------------------------------------------
        if all(v2 > 9):
            self.N    = {x: Tab1[x] for x in bds}
            self.H    = {x: Tab2[x] for x in bds}
            self.NH   = {x: [Tab1[x], Tab2[x]] for x in bds}
            self.NHHz = {x: [Tab1[x] * 60, Tab2[x] * 600] for x in bds}
            for x in bds:
                if x[0] in ('G', 'A'):
                    self.NHnu[x] = {'N1': Tab1[x], 'H1': Tab2[x]}
                if x[0] in ('U', 'C'):
                    self.NHnu[x] = {'N3': Tab1[x], 'H3': Tab2[x]}

        # -------- 13C / 1H (aromatic C8/H8, C6/H6) ---------------------------
        elif all(v2 < 8.6) and all(v2 > 6.4) and all(v1 < 148) and all(v1 > 134):
            self.C    = {x: Tab1[x] for x in bds}
            self.H    = {x: Tab2[x] for x in bds}
            self.CH   = {x: [Tab1[x], Tab2[x]] for x in bds}
            self.CHHz = {x: [Tab1[x] * 150, Tab2[x] * 600] for x in bds}
            for x in bds:
                if x[0] in ('G', 'A'):
                    self.CHnu[x] = {'C8': Tab1[x], 'H8': Tab2[x]}
                if x[0] in ('U', 'C'):
                    self.CHnu[x] = {'C6': Tab1[x], 'H6': Tab2[x]}

        # -------- 13C / 1H (sugar C1'/H1') ----------------------------------
        elif all(v2 < 6.4) and all(v2 > 4.8) and all(v1 < 98) and all(v1 > 94):
            self.C    = {x: Tab1[x] for x in bds}
            self.H    = {x: Tab2[x] for x in bds}
            self.CH   = {x: [Tab1[x], Tab2[x]] for x in bds}
            self.CHHz = {x: [Tab1[x] * 150, Tab2[x] * 600] for x in bds}
            for x in bds:
                self.CHnu[x] = {"C1'": Tab1[x], "H1'": Tab2[x]}

        # -------- 13C / 1H (aromatic C2/H2 of adenine) ----------------------
        elif all(v2 < 8.4) and all(v2 > 6.0) and all(v1 < 158) and all(v1 > 148):
            self.C    = {x: Tab1[x] for x in bds}
            self.H    = {x: Tab2[x] for x in bds}
            self.CH   = {x: [Tab1[x], Tab2[x]] for x in bds}
            self.CHHz = {x: [Tab1[x] * 150, Tab2[x] * 600] for x in bds}
            for x in bds:
                self.CHnu[x] = {"C2": Tab1[x], "H2": Tab2[x]}
