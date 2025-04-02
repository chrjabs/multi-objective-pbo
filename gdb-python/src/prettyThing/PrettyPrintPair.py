import gdb.printing
import gdb.types

# taken from /usr/share/gcc-11/python/libstdcxx/v6/printers.py
Iterator = object


class PairPrinter:

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self):
        return ('{%s, %s}'
                % (self.val['first'], self.val['second']))

    def display_hint(self):
        return 'array'
