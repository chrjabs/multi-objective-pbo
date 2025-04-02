class LboolPrinter:
    "Print a Lit"

    def __init__(self, val):
        self.val = val
        if self.val['value'] == 0:
            self.pol = True
        else:
            self.pol = False

    def to_string(self):
        if self.pol:
            return str(1)
        else:
            return str(0)
