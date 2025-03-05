###############################################################################
# LLDB Script to improve introspection of array types when debugging software
# using Enoki. Copy this file to "~/.lldb" (creating the directory, if not
# present) and then apppend the following line to the file "~/.lldbinit"
# (again, creating it if, not already present):
###############################################################################
# command script import ~/.lldb/enoki_lldb.py
###############################################################################

import sys
import lldb

simple_types = {
    'bool',
    'char', 'unsigned char',
    'short', 'unsigned short',
    'int', 'unsigned int',
    'long', 'unsigned long',
    'long long', 'unsigned long long',
    'float', 'double'
}


class StaticArraySynthProvider:
    def __init__(self, instance, internal_dict):
        self.instance = instance

    def update(self):
        itype = self.instance.GetType().GetCanonicalType().GetUnqualifiedType()
        itype_name = itype.name

        # Extract derived type
        if itype_name.startswith('enoki::StaticArrayImpl'):
            itype = itype.GetTemplateArgumentType(3)
            itype_name = itype.name

        # Determine the size
        self.size = int(itype_name[itype_name.rfind(',')+1:itype_name.rfind('>')])

        self.is_mask = 'Mask' in itype_name

        data = self.instance.GetChildMemberWithName('m_data')
        if data:
            self.data_type = data.GetType().GetTemplateArgumentType(0)
        else:
            self.data_type = itype.GetTemplateArgumentType(0)

        self.type_size = self.data_type.GetByteSize()
        self.is_simple = self.data_type.name in simple_types
        self.kmask = self.instance.GetChildMemberWithName('k')

    def has_children(self):
        return not self.is_simple and self.size > 0

    def num_children(self):
        return 0 if self.is_simple else self.size

    def get_child_index(self, name):
        try:
            return int(name)
        except Exception:
            return -1

    def get_child_at_index(self, index):
        if index < 0 or index >= self.size:
            return None
        return self.instance.CreateChildAtOffset(
            str(index), index * self.type_size, self.data_type)

    def get_summary(self):
        if self.is_simple:
            if not self.is_mask:
                result = [str(self.get_child_at_index(i).value) for i in range(self.size)]
            else:
                if self.kmask:
                    # AVX512 mask register
                    result = list(reversed(format(int(self.kmask.unsigned), '0%ib' % self.size)))
                else:
                    result = [None] * self.size
                    for i in range(self.size):
                        value = self.get_child_at_index(i).value
                        result[i] = '0' if (value == '0' or value == 'false') else '1'
            return '[' + ', '.join(result) + ']'
        else:
            return ''


class DynamicArraySummaryProvider:
    def __init__(self, instance, internal_dict):
        self.instance = instance

    def update(self):
        self.size = self.instance.GetChildMemberWithName('m_size').unsigned
        self.packet_count = self.instance.GetChildMemberWithName('m_packets_allocated').unsigned
        self.packet_type = self.instance.GetType().GetCanonicalType().\
            GetUnqualifiedType().GetTemplateArgumentType(0)
        self.packet_size = self.packet_type.GetByteSize()
        self.ptr = self.instance.GetChildMemberWithName('m_packets').GetData()
        error = lldb.SBError()
        self.ptr = self.ptr.GetUnsignedInt64(offset=0, error=error) if self.ptr.GetByteSize() == 8 \
              else self.ptr.GetUnsignedInt32(offset=0, error=error)
        self.limit = 20

    def has_children(self):
        return False

    def num_children(self):
        return 0

    def get_child_index(self, name):
        return None

    def get_child_at_index(self, index):
        return None

    def get_summary(self):
        values = []
        for i in range(self.packet_count):
            value = str(self.instance.CreateValueFromAddress(str(i),
                self.ptr + i*self.packet_size, self.packet_type))
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


def attach(enoki_category, synth_class, type_name, summary=True, synth=True):
    if summary:
        def summary_func(instance, internal_dict):
            synth = synth_class(instance.GetNonSyntheticValue(), internal_dict)
            synth.update()
            return synth.get_summary()

        summary_func.__name__ = synth_class.__name__ + 'SummaryWrapper'
        setattr(sys.modules[__name__], summary_func.__name__, summary_func)

        summary = lldb.SBTypeSummary.CreateWithFunctionName(__name__ + '.' + summary_func.__name__)
        summary.SetOptions(lldb.eTypeOptionCascade)
        enoki_category.AddTypeSummary(lldb.SBTypeNameSpecifier(type_name, True), summary)

    if synth:
        synth = lldb.SBTypeSynthetic.CreateWithClassName(__name__ + '.' + synth_class.__name__)
        synth.SetOptions(lldb.eTypeOptionCascade)
        enoki_category.AddTypeSynthetic(lldb.SBTypeNameSpecifier(type_name, True), synth)


def __lldb_init_module(debugger, internal_dict):
    enoki_category = debugger.CreateCategory('enoki')
    enoki_category.SetEnabled(True)

    # Static Enoki arrays
    regexp_1 = r'enoki::(Array|Packet|Complex|Matrix|' \
        'Quaternion|StaticArrayImpl)(Mask)?<.+>'

    # Mitsuba 2 is one of the main users of Enoki. For convenience, also
    # declare its custom array types here
    regexp_2 = r'mitsuba::(Vector|Point|Normal|Spectrum|Color)<.+>'

    regexp_combined = r'^(%s)|(%s)$' % (regexp_1, regexp_2)
    attach(enoki_category, StaticArraySynthProvider, regexp_combined)

    # Dynamic Enoki arrays
    attach(enoki_category, DynamicArraySummaryProvider,
           r"^enoki::DynamicArray(Impl)?<.+>$")
