import gdb.printing
import gdb.types

# taken from /usr/share/gcc-11/python/libstdcxx/v6/printers.py
Iterator = object


class YPointPrinter:
    "Print a yPoint"

    class _iterator(Iterator):
        def __init__(self, start, finish):
            self.item = start
            self.finish = finish
            self.count = 0

        def __iter__(self):
            return self

        def __next__(self):
            count = self.count
            self.count = self.count + 1
            if self.item == self.finish:
                raise StopIteration
            elt = self.item.dereference()
            self.item = self.item + 1
            return ('f[%d] -> %d' % (count, elt))

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self):
        start = self.val['_M_impl']['_M_start']
        finish = self.val['_M_impl']['_M_finish']
        it = self._iterator(start, finish)
        children = ', '.join([string for string in it])
        return ('{t->%s, i->{l->%d}, d->{%s}}'
                % (self.val.type.tag, int(finish - start), children))

    def display_hint(self):
        return
