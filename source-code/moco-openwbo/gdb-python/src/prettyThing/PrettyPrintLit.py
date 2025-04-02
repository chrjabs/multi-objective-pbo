import gdb.printing
import gdb.types


class LitPrinter:
    "Print a Lit"

    def __init__(self, val):
        pass
        self.val = val

    def to_string(self):
        def sign(val):
            return val['x'] & 1

        def var(val):
            return val['x'] >> 1

        if sign(self.val):
            # +1, so that inner and outer representations match
            return '~x%d' % (var(self.val) + 1)
        else:
            return 'x%d' % (var(self.val) + 1)
