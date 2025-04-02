import gdb.printing
import gdb.types

# taken from /usr/share/gcc-11/python/libstdcxx/v6/printers.py
Iterator = object


class PBObjFunctionPrinter:
    "Print a PBObjFunction"

    class _iterator(Iterator):
        def __init__(self, first_lits, first_coeffs, length):
            self.lit_item = first_lits
            self.coeff_item = first_coeffs
            self.length = length
            self.min = int(self.coeff_item.dereference())
            self.max = self.min
            self.total = 0
            self.count = 0

        def __iter__(self):
            return self

        def __next__(self):
            if self.count == self.length:
                raise StopIteration
            elt_lit = self.lit_item.dereference()
            elt_coeff = self.coeff_item.dereference()
            coeff = int(elt_coeff)
            self.total += coeff
            if self.min > coeff:
                self.min = coeff
            if self.max < coeff:
                self.max = coeff
            self.count = self.count + 1
            self.lit_item = self.lit_item + 1
            self.coeff_item = self.coeff_item + 1
            return ('%+d %s' %
                    (elt_coeff, elt_lit))

    def __init__(self, val: gdb.Value):
        self.val = val
        self.length = self.val['_lits']['sz']

    def to_string(self):
        it = self._iterator(self.val['_lits']['data'],
                            self.val['_coeffs']['data'],
                            self.val['_lits']['sz'])
        children = ' '.join([string for string in it])
        return '{t->%s, i->%s, d->{%s}}' % (
            self.val.type,
            '{l->%d, max->%d, min->%d, total->%d}' % (
                self.length,
                it.max,
                it.min,
                it.total),
            children)
