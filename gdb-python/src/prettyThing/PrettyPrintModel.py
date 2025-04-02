import gdb.printing
import gdb.types

# taken from /usr/share/gcc-11/python/libstdcxx/v6/printers.py
Iterator = object


class ModelPrinter:
    "Print a PBObjFunction"

    class _iterator(Iterator):
        def __init__(self, first, last):
            self.item = first
            self.last = last
            self.count = 0

        def __iter__(self):
            return self

        def __next__(self):
            if self.item == self.last:
                raise StopIteration
            # elt will be printed according to the lbool printer
            elt = self.item.dereference()
            # +1, so that external x1 is printed like x1
            var = self.count + 1
            self.count = self.count + 1
            self.item = self.item + 1
            return 'x%d -> %s' % (var, elt)

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self):
        start = self.val['_M_impl']['_M_start']
        end = self.val['_M_impl']['_M_finish']
        it = self._iterator(start,
                            end)

        return '{t->%s, i->%s, d->{%s}}' % (
            self.val.type.name,
            '{l->%d}' % (
                int(end - start)),
            ', '.join([string for string in it]))
