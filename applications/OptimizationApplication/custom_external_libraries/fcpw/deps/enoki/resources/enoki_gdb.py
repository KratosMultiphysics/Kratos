###############################################################################
# GDB Script to improve introspection of array types when debugging software
# using Enoki. Copy this file to "~/.gdb" (creating the directory, if not
# present) and then apppend the following line to the file "~/.gdbinit"
# (again, creating it if, not already present):
###############################################################################
# set print pretty
# source ~/.gdb/enoki_gdb.py
###############################################################################

import gdb

simple_types = {
    'bool',
    'char', 'unsigned char',
    'short', 'unsigned short',
    'int', 'unsigned int',
    'long', 'unsigned long',
    'long long', 'unsigned long long',
    'float', 'double'
}


class EnokiIterator:
    def __init__(self, instance, size):
        self.instance = instance
        self.size = size
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= self.size:
            raise StopIteration
        result = ('[%i]' % self.index, self.instance[self.index])
        self.index += 1
        return result

    def next(self):
        return self.__next__()


class EnokiStaticArrayPrinter:
    def __init__(self, instance):
        self.instance = instance
        itype = self.instance.type.strip_typedefs()

        # Extract derived type
        if 'StaticArrayImpl' in str(itype):
            itype = itype.template_argument(3)

        try:
            data = self.instance['m_data']['_M_elems']
            self.data_type = data.type.strip_typedefs().target()
        except Exception:
            self.data_type = itype.template_argument(0)

        # Determine the size and data type
        self.size = int(str(itype.template_argument(1)))
        self.is_simple = str(self.data_type) in simple_types
        self.type_size = self.data_type.sizeof
        self.is_mask = 'Mask' in str(itype)

        try:
            _ = instance['k']
            self.kmask = True
        except Exception:
            self.kmask = False

    def entry(self, i):
        if i < 0 or i >= self.size:
            return None
        addr = int(self.instance.address) + self.type_size * i
        cmd = '*((%s *) 0x%x)' % (str(self.data_type), addr)
        return str(gdb.parse_and_eval(cmd))

    def children(self):
        if self.is_simple:
            return []
        else:
            return EnokiIterator(self.instance['m_data']['_M_elems'], self.size)

    def to_string(self):
        if self.is_simple:
            if not self.is_mask:
                result = [self.entry(i) for i in range(self.size)]
            else:
                if self.kmask:
                    # AVX512 mask register
                    result = list(reversed(format(int(self.instance['k']), '0%ib' % self.size)))
                else:
                    result = [None] * self.size
                    for i in range(self.size):
                        value = self.entry(i)
                        result[i] = '0' if (value == '0' or value == 'false') else '1'
            return '[' + ', '.join(result) + ']'
        else:
            return ''


class EnokiDynamicArrayPrinter:
    def __init__(self, instance):
        self.instance = instance
        itype = self.instance.type.strip_typedefs()
        self.size = int(str(self.instance['m_size']))
        self.packet_count = int(str(self.instance['m_packets_allocated']))
        self.packet_type = itype.template_argument(0)
        self.packet_size = self.packet_type.sizeof
        self.data = int(str(instance['m_packets']['_M_t']['_M_t']['_M_head_impl']), 0)
        self.limit = 20

    def to_string(self):
        values = []
        for i in range(self.packet_count):
            addr = int(self.data) + self.packet_size * i
            cmd = '*((%s *) 0x%x)' % (str(self.packet_type), addr)
            value = str(gdb.parse_and_eval(cmd))
            assert value[-1] == ']'
            values += value[value.rfind('[')+1:-1].split(', ')
            if len(values) > self.size:
                values = values[0:self.size]
                break
            if len(values) > self.limit:
                break
        if len(values) > self.limit:
            values = values[0:self.limit]
            values.append(".. %i skipped .." % (self.size - self.limit))
        return '[' + ', '.join(values) + ']'

# Static Enoki arrays
regexp_1 = r'(enoki::)?(Array|Packet|Complex|Matrix|' \
    'Quaternion|StaticArrayImpl)(Mask)?<.+>'

# Mitsuba 2 is one of the main users of Enoki. For convenience, also
# declare its custom array types here
regexp_2 = r'(mitsuba::)?(Vector|Point|Normal|Spectrum|Color)<.+>'

regexp_combined = r'^(%s)|(%s)$' % (regexp_1, regexp_2)

p = gdb.printing.RegexpCollectionPrettyPrinter("enoki")
p.add_printer("static", regexp_combined, EnokiStaticArrayPrinter)
p.add_printer("dynamic", r'^(enoki::)?DynamicArray(Impl)?<.+>$', EnokiDynamicArrayPrinter)

o = gdb.current_objfile()
gdb.printing.register_pretty_printer(o, p)
