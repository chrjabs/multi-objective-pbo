import gdb.printing
import gdb.types

# taken from /usr/share/gcc-11/python/libstdcxx/v6/printers.py
Iterator = object


class StdVecPrinter:
    "Print a std::vector"

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
            return ('[%d]' % count, elt)

    def __init__(self, val: gdb.Value):
        self.val = val
        self.inner_type = self.val.type.template_argument(0)
        self.length = self.val['sz']

    def children(self):
        start = self.val['data']
        finish = start + self.length
        return self._iterator(start, finish)

    def to_string(self):

        capacity = self.val['cap']
        return ('{t->vecOf%s, i->{l->%d, c -> %d}}'
                % (self.inner_type, self.length, int(capacity)))

    def display_hint(self):
        return 'array'
