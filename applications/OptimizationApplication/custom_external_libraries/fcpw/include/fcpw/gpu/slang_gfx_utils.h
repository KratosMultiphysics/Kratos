#pragma once

#include "slang-gfx.h"

namespace fcpw {

using namespace gfx;

class GPUContext {
public:
    IDevice::Desc deviceDesc = {};
    ComPtr<IDevice> device;
    ITransientResourceHeap::Desc transientHeapDesc = {};
    ComPtr<ITransientResourceHeap> transientHeap;
    ICommandQueue::Desc queueDesc = { ICommandQueue::QueueType::Graphics };
    ComPtr<ICommandQueue> queue;

    void initDevice(slang::PreprocessorMacroDesc& macro, int nMacros) {
        deviceDesc.slang.preprocessorMacros = &macro;
        deviceDesc.slang.preprocessorMacroCount = nMacros;
        SlangResult createDeviceResult = gfxCreateDevice(&deviceDesc, device.writeRef());
        if (createDeviceResult != SLANG_OK) {
            std::cout << "failed to create device" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::cout << "device: " << device->getDeviceInfo().apiName << std::endl;
    }

    void initTransientResources() {
        transientHeapDesc.constantBufferSize = 4096;
        SlangResult createTransientHeapResult = device->createTransientResourceHeap(
            transientHeapDesc, transientHeap.writeRef());
        if (createTransientHeapResult != SLANG_OK) {
            std::cout << "failed to create transient resource heap" << std::endl;
            exit(EXIT_FAILURE);
        }

        queue = device->createCommandQueue(queueDesc);
    }
};

class Shader {
public:
    ComPtr<IShaderProgram> program;
    slang::ProgramLayout* reflection = nullptr;
    ComputePipelineStateDesc pipelineDesc = {};
    ComPtr<IPipelineState> pipelineState;

    Slang::Result loadComputeProgram(ComPtr<IDevice>& device,
                                     const char* shaderModuleName,
                                     const char* entryPointName) {
        // load shader module
        ComPtr<slang::ISession> slangSession;
        SLANG_RETURN_ON_FAIL(device->getSlangSession(slangSession.writeRef()));
        ComPtr<slang::IBlob> diagnosticsBlob;
        slang::IModule* module = slangSession->loadModule(shaderModuleName, diagnosticsBlob.writeRef());
        diagnoseIfNeeded(diagnosticsBlob);
        if (!module) return SLANG_FAIL;

        // set entry point and create composite program
        ComPtr<slang::IEntryPoint> computeEntryPoint;
        SLANG_RETURN_ON_FAIL(module->findEntryPointByName(entryPointName, computeEntryPoint.writeRef()));

        slang::IComponentType* componentTypes[] = { module, computeEntryPoint };
        ComPtr<slang::IComponentType> composedProgram;
        SlangResult result = slangSession->createCompositeComponentType(componentTypes, 2,
                                                                        composedProgram.writeRef(),
                                                                        diagnosticsBlob.writeRef());
        diagnoseIfNeeded(diagnosticsBlob);
        SLANG_RETURN_ON_FAIL(result);

        // link program and get reflection
        ComPtr<slang::IComponentType> linkedProgram;
        result = composedProgram->link(linkedProgram.writeRef(), diagnosticsBlob.writeRef());
        diagnoseIfNeeded(diagnosticsBlob);
        SLANG_RETURN_ON_FAIL(result);

        composedProgram = linkedProgram;
        reflection = composedProgram->getLayout();

        // create shader pipeline state
        IShaderProgram::Desc programDesc = {};
        programDesc.slangGlobalScope = composedProgram.get();

        program = device->createProgram(programDesc);
        pipelineDesc.program = program.get();

        SLANG_RETURN_ON_FAIL(device->createComputePipelineState(pipelineDesc, pipelineState.writeRef()));
        return SLANG_OK;
    }

private:
    void diagnoseIfNeeded(slang::IBlob* diagnosticsBlob) {
        if (diagnosticsBlob != nullptr) {
            std::cerr << (const char*)diagnosticsBlob->getBufferPointer() << std::endl;
        }
    }
};

class GPUBuffer {
public:
    IBufferResource::Desc desc;
    ComPtr<IBufferResource> buffer;
    ComPtr<IResourceView> view;

    template<typename T>
    Slang::Result create(ComPtr<IDevice>& device, bool unorderedAccess,
                         T* initialData, size_t elementCount) {
        T *data = nullptr;
        if (elementCount == 0) {
            elementCount = 1; // Slang requires buffers to be non-empty

        } else {
            data = initialData;
        }

        desc.sizeInBytes = elementCount * sizeof(T);
        desc.format = Format::Unknown;
        desc.elementSize = sizeof(T);
        desc.defaultState = unorderedAccess ? ResourceState::UnorderedAccess :
                                              ResourceState::ShaderResource;
        desc.memoryType = MemoryType::DeviceLocal;
        desc.allowedStates = ResourceStateSet(ResourceState::ShaderResource,
                                              ResourceState::CopyDestination,
                                              ResourceState::CopySource);
        if (unorderedAccess) desc.allowedStates.add(ResourceState::UnorderedAccess);
        SLANG_RETURN_ON_FAIL(device->createBufferResource(desc, (void *)data, buffer.writeRef()));

        IResourceView::Desc viewDesc = {};
        viewDesc.type = unorderedAccess ? IResourceView::Type::UnorderedAccess :
                                          IResourceView::Type::ShaderResource;
        viewDesc.format = Format::Unknown;
        SLANG_RETURN_ON_FAIL(device->createBufferView(buffer, nullptr, viewDesc, view.writeRef()));

        return SLANG_OK;
    }

    template<typename T>
    Slang::Result read(ComPtr<IDevice>& device, size_t elementCount, std::vector<T>& result) const {
        ComPtr<ISlangBlob> resultBlob;
        size_t expectedBufferSize = elementCount * sizeof(T);
        SLANG_RETURN_ON_FAIL(device->readBufferResource(buffer, 0, expectedBufferSize, resultBlob.writeRef()));
        if (resultBlob->getBufferSize() != expectedBufferSize) {
            std::cout << "incorrect GPU buffer size on read" << std::endl;
            return SLANG_FAIL;
        }

        result.clear();
        auto resultPtr = (T *)resultBlob->getBufferPointer();
        result.assign(resultPtr, resultPtr + elementCount);

        return SLANG_OK;
    }
};

/// Represents a "pointer" to the storage for a shader parameter of a (dynamically) known type.
///
/// A `ShaderCursor` serves as a pointer-like type for things stored inside a `ShaderObject`.
///
/// A cursor that points to the entire content of a shader object can be formed as
/// `ShaderCursor(someObject)`. A cursor pointing to a structure field or array element can be
/// formed from another cursor using `getField` or `getElement` respectively.
///
/// Given a cursor pointing to a value of some "primitive" type, we can set or get the value
/// using operations like `setResource`, `getResource`, etc.
///
/// Because type information for shader parameters is being reflected dynamically, all type
/// checking for shader cursors occurs at runtime, and errors may occur when attempting to
/// set a parameter using a value of an inappropriate type. As much as possible, `ShaderCursor`
/// attempts to protect against these cases and return an error `Result` or an invalid
/// cursor, rather than allowing operations to proceed with incorrect types.
///
struct ShaderCursor
{
    IShaderObject* m_baseObject = nullptr;
    slang::TypeLayoutReflection* m_typeLayout = nullptr;
    ShaderObjectContainerType m_containerType = ShaderObjectContainerType::None;
    ShaderOffset m_offset;

    /// Get the type (layout) of the value being pointed at by the cursor
    slang::TypeLayoutReflection* getTypeLayout() const { return m_typeLayout; }

    /// Is this cursor valid (that is, does it seem to point to an actual location)?
    ///
    /// This check is equivalent to checking whether a pointer is null, so it is
    /// a very weak sense of "valid." In particular, it is possible to form a
    /// `ShaderCursor` for which `isValid()` is true, but attempting to get or
    /// set the value would be an error (like dereferencing a garbage pointer).
    ///
    bool isValid() const { return m_baseObject != nullptr; }

    Result getDereferenced(ShaderCursor& outCursor) const;

    ShaderCursor getDereferenced() const
    {
        ShaderCursor result;
        getDereferenced(result);
        return result;
    }

    /// Form a cursor pointing to the field with the given `name` within the value this cursor
    /// points at.
    ///
    /// If the operation succeeds, then the field cursor is written to `outCursor`.
    Result getField(const char* nameBegin, const char* nameEnd, ShaderCursor& outCursor) const;

    ShaderCursor getField(const char* name) const
    {
        ShaderCursor cursor;
        getField(name, nullptr, cursor);
        return cursor;
    }

    /// Some resources such as RWStructuredBuffer, AppendStructuredBuffer and
    /// ConsumeStructuredBuffer need to have their counter explicitly bound on
    /// APIs other than DirectX, this will return a valid ShaderCursor pointing
    /// to that resource if that is the case.
    /// Otherwise, this returns an invalid cursor.
    ShaderCursor getExplicitCounter() const;

    ShaderCursor getElement(GfxIndex index) const;

    static Result followPath(const char* path, ShaderCursor& ioCursor);

    ShaderCursor getPath(const char* path) const
    {
        ShaderCursor result(*this);
        followPath(path, result);
        return result;
    }

    ShaderCursor() {}

    ShaderCursor(IShaderObject* object)
        : m_baseObject(object)
        , m_typeLayout(object->getElementTypeLayout())
        , m_containerType(object->getContainerType())
    {}

    SlangResult setData(void const* data, Size size) const
    {
        return m_baseObject->setData(m_offset, data, size);
    }

    template <typename T>
    SlangResult setData(T const& data) const
    {
        return setData(&data, sizeof(data));
    }

    SlangResult setObject(IShaderObject* object) const
    {
        return m_baseObject->setObject(m_offset, object);
    }

    SlangResult setSpecializationArgs(const slang::SpecializationArg* args, GfxCount count) const
    {
        return m_baseObject->setSpecializationArgs(m_offset, args, count);
    }

    SlangResult setResource(IResourceView* resourceView) const
    {
        return m_baseObject->setResource(m_offset, resourceView);
    }

    SlangResult setSampler(ISamplerState* sampler) const
    {
        return m_baseObject->setSampler(m_offset, sampler);
    }

    SlangResult setCombinedTextureSampler(IResourceView* textureView, ISamplerState* sampler) const
    {
        return m_baseObject->setCombinedTextureSampler(m_offset, textureView, sampler);
    }

        /// Produce a cursor to the field with the given `name`.
        ///
        /// This is a convenience wrapper around `getField()`.
    ShaderCursor operator[](const char* name) const
    {
        return getField(name);
    }

        /// Produce a cursor to the element or field with the given `index`.
        ///
        /// This is a convenience wrapper around `getElement()`.
    ShaderCursor operator[](int64_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](uint64_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](int32_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](uint32_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](int16_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](uint16_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](int8_t index) const { return getElement((GfxIndex)index); }
    ShaderCursor operator[](uint8_t index) const { return getElement((GfxIndex)index); }
};

Result ShaderCursor::getDereferenced(ShaderCursor& outCursor) const
{
    switch (m_typeLayout->getKind())
    {
    default:
        return SLANG_E_INVALID_ARG;

    case slang::TypeReflection::Kind::ConstantBuffer:
    case slang::TypeReflection::Kind::ParameterBlock:
        {
            auto subObject = m_baseObject->getObject(m_offset);
            outCursor = ShaderCursor(subObject);
            return SLANG_OK;
        }
    }
}

ShaderCursor ShaderCursor::getExplicitCounter() const
{
    // Similar to getField below

    // The alternative to handling this here would be to augment IResourceView
    // with a `getCounterResourceView()`, and set that also in `setResource`
    if(const auto counterVarLayout = m_typeLayout->getExplicitCounter())
    {
        ShaderCursor counterCursor;

        // The counter cursor will point into the same parent object.
        counterCursor.m_baseObject = m_baseObject;

        // The type being pointed to is the type of the field.
        counterCursor.m_typeLayout = counterVarLayout->getTypeLayout();

        // The byte offset is the current offset plus the relative offset of the counter.
        // The offset in binding ranges is computed similarly.
        counterCursor.m_offset.uniformOffset
            = m_offset.uniformOffset + SlangInt(counterVarLayout->getOffset());
        counterCursor.m_offset.bindingRangeIndex
            = m_offset.bindingRangeIndex + GfxIndex(m_typeLayout->getExplicitCounterBindingRangeOffset());

        // The index of the counter within any binding ranges will be the same
        // as the index computed for the parent structure.
        //
        // Note: this case would arise for an array of structured buffers
        //
        //      AppendStructuredBuffer g[4];
        //
        // In this scenario, `g` holds two binding ranges:
        //
        // * Range #0 comprises 4 element buffers, representing `g[...].elements`
        // * Range #1 comprises 4 counter buffers, representing `g[...].counter`
        //
        // A cursor for `g[2]` would have a `bindingRangeIndex` of zero but
        // a `bindingArrayIndex` of 2, indicating that we could end up
        // referencing either range, but no matter what we know the index
        // is 2. Thus when we form a cursor for `g[2].counter` we want to
        // apply the binding range offset to get a `bindingRangeIndex` of
        // 1, while the `bindingArrayIndex` is unmodified.
        //
        // The result is that `g[2].counter` is stored in range #1 at array index 2.
        //
        counterCursor.m_offset.bindingArrayIndex = m_offset.bindingArrayIndex;

        return counterCursor;
    }
    // Otherwise, return an invalid cursor
    return ShaderCursor{};
}

Result ShaderCursor::getField(const char* name, const char* nameEnd, ShaderCursor& outCursor) const
{
    // If this cursor is invalid, then can't possible fetch a field.
    //
    if (!isValid())
        return SLANG_E_INVALID_ARG;

    // If the cursor is valid, we want to consider the type of data
    // it is referencing.
    //
    switch (m_typeLayout->getKind())
    {
        // The easy/expected case is when the value has a structure type.
        //
    case slang::TypeReflection::Kind::Struct:
        {
            // We start by looking up the index of a field matching `name`.
            //
            // If there is no such field, we have an error.
            //
            SlangInt fieldIndex = m_typeLayout->findFieldIndexByName(name, nameEnd);
            if (fieldIndex == -1)
                break;

            // Once we know the index of the field being referenced,
            // we create a cursor to point at the field, based on
            // the offset information already in this cursor, plus
            // offsets derived from the field's layout.
            //
            slang::VariableLayoutReflection* fieldLayout =
                m_typeLayout->getFieldByIndex((unsigned int)fieldIndex);
            ShaderCursor fieldCursor;

            // The field cursorwill point into the same parent object.
            //
            fieldCursor.m_baseObject = m_baseObject;

            // The type being pointed to is the tyep of the field.
            //
            fieldCursor.m_typeLayout = fieldLayout->getTypeLayout();

            // The byte offset is the current offset plus the relative offset of the field.
            // The offset in binding ranges is computed similarly.
            //
            fieldCursor.m_offset.uniformOffset = m_offset.uniformOffset + fieldLayout->getOffset();
            fieldCursor.m_offset.bindingRangeIndex =
                m_offset.bindingRangeIndex + (GfxIndex)m_typeLayout->getFieldBindingRangeOffset(fieldIndex);

            // The index of the field within any binding ranges will be the same
            // as the index computed for the parent structure.
            //
            // Note: this case would arise for an array of structures with texture-type
            // fields. Suppose we have:
            //
            //      struct S { Texture2D t; Texture2D u; }
            //      S g[4];
            //
            // In this scenario, `g` holds two binding ranges:
            //
            // * Range #0 comprises 4 textures, representing `g[...].t`
            // * Range #1 comprises 4 textures, representing `g[...].u`
            //
            // A cursor for `g[2]` would have a `bindingRangeIndex` of zero but
            // a `bindingArrayIndex` of 2, iindicating that we could end up
            // referencing either range, but no matter what we know the index
            // is 2. Thus when we form a cursor for `g[2].u` we want to
            // apply the binding range offset to get a `bindingRangeIndex` of
            // 1, while the `bindingArrayIndex` is unmodified.
            //
            // The result is that `g[2].u` is stored in range #1 at array index 2.
            //
            fieldCursor.m_offset.bindingArrayIndex = m_offset.bindingArrayIndex;

            outCursor = fieldCursor;
            return SLANG_OK;
        }
        break;

    // In some cases the user might be trying to acess a field by name
    // from a cursor that references a constant buffer or parameter block,
    // and in these cases we want the access to Just Work.
    //
    case slang::TypeReflection::Kind::ConstantBuffer:
    case slang::TypeReflection::Kind::ParameterBlock:
        {
            // We basically need to "dereference" the current cursor
            // to go from a pointer to a constant buffer to a pointer
            // to the *contents* of the constant buffer.
            //
            ShaderCursor d = getDereferenced();
            return d.getField(name, nameEnd, outCursor);
        }
        break;
    }

    // If a cursor is pointing at a root shader object (created for a
    // program), then we will also iterate over the entry point shader
    // objects attached to it and look for a matching parameter name
    // on them.
    //
    // This is a bit of "do what I mean" logic and could potentially
    // lead to problems if there could be multiple entry points with
    // the same parameter name.
    //
    // TODO: figure out whether we should support this long-term.
    //
    auto entryPointCount = (GfxIndex) m_baseObject->getEntryPointCount();
    for( GfxIndex e = 0; e < entryPointCount; ++e )
    {
        ComPtr<IShaderObject> entryPoint;
        m_baseObject->getEntryPoint(e, entryPoint.writeRef());

        ShaderCursor entryPointCursor(entryPoint);

        auto result = entryPointCursor.getField(name, nameEnd, outCursor);
        if(SLANG_SUCCEEDED(result))
            return result;
    }

    return SLANG_E_INVALID_ARG;
}

ShaderCursor ShaderCursor::getElement(GfxIndex index) const
{
    if (m_containerType != ShaderObjectContainerType::None)
    {
        ShaderCursor elementCursor;
        elementCursor.m_baseObject = m_baseObject;
        elementCursor.m_typeLayout = m_typeLayout->getElementTypeLayout();
        elementCursor.m_containerType = m_containerType;
        elementCursor.m_offset.uniformOffset = index * m_typeLayout->getStride();
        elementCursor.m_offset.bindingRangeIndex = 0;
        elementCursor.m_offset.bindingArrayIndex = index;
        return elementCursor;
    }

    switch( m_typeLayout->getKind() )
    {
    case slang::TypeReflection::Kind::Array:
        {
            ShaderCursor elementCursor;
            elementCursor.m_baseObject = m_baseObject;
            elementCursor.m_typeLayout = m_typeLayout->getElementTypeLayout();
            elementCursor.m_offset.uniformOffset =
                m_offset.uniformOffset +
                index * m_typeLayout->getElementStride(SLANG_PARAMETER_CATEGORY_UNIFORM);
            elementCursor.m_offset.bindingRangeIndex = m_offset.bindingRangeIndex;
            elementCursor.m_offset.bindingArrayIndex =
                m_offset.bindingArrayIndex * (GfxCount)m_typeLayout->getElementCount() + index;
            return elementCursor;
        }
        break;

    case slang::TypeReflection::Kind::Struct:
        {
            // The logic here is similar to `getField()` except that we don't
            // need to look up the field index based on a name first.
            //
            auto fieldIndex = index;
            slang::VariableLayoutReflection* fieldLayout =
                m_typeLayout->getFieldByIndex((unsigned int)fieldIndex);
            if(!fieldLayout)
                return ShaderCursor();

            ShaderCursor fieldCursor;
            fieldCursor.m_baseObject = m_baseObject;
            fieldCursor.m_typeLayout = fieldLayout->getTypeLayout();
            fieldCursor.m_offset.uniformOffset = m_offset.uniformOffset + fieldLayout->getOffset();
            fieldCursor.m_offset.bindingRangeIndex =
                m_offset.bindingRangeIndex + (GfxIndex)m_typeLayout->getFieldBindingRangeOffset(fieldIndex);
            fieldCursor.m_offset.bindingArrayIndex = m_offset.bindingArrayIndex;

            return fieldCursor;
        }
        break;

    case slang::TypeReflection::Kind::Vector:
    case slang::TypeReflection::Kind::Matrix:
        {
            ShaderCursor fieldCursor;
            fieldCursor.m_baseObject = m_baseObject;
            fieldCursor.m_typeLayout = m_typeLayout->getElementTypeLayout();
            fieldCursor.m_offset.uniformOffset = m_offset.uniformOffset + m_typeLayout->getElementStride(SLANG_PARAMETER_CATEGORY_UNIFORM) * index;
            fieldCursor.m_offset.bindingRangeIndex = m_offset.bindingRangeIndex;
            fieldCursor.m_offset.bindingArrayIndex = m_offset.bindingArrayIndex;
            return fieldCursor;
        }
        break;
    }

    return ShaderCursor();
}


static int _peek(const char* slice)
{
    const char* b = slice;
    if (!b || !*b)
        return -1;
    return *b;
}

static int _get(const char*& slice)
{
    const char* b = slice;
    if (!b || !*b)
        return -1;
    auto result = *b++;
    slice = b;
    return result;
}

Result ShaderCursor::followPath(const char* path, ShaderCursor& ioCursor)
{
    ShaderCursor cursor = ioCursor;

    enum
    {
        ALLOW_NAME = 0x1,
        ALLOW_SUBSCRIPT = 0x2,
        ALLOW_DOT = 0x4,
    };
    int state = ALLOW_NAME | ALLOW_SUBSCRIPT;

    const char* rest = path;
    for (;;)
    {
        int c = _peek(rest);

        if (c == -1)
            break;
        else if (c == '.')
        {
            if (!(state & ALLOW_DOT))
                return SLANG_E_INVALID_ARG;

            _get(rest);
            state = ALLOW_NAME;
            continue;
        }
        else if (c == '[')
        {
            if (!(state & ALLOW_SUBSCRIPT))
                return SLANG_E_INVALID_ARG;

            _get(rest);
            GfxCount index = 0;
            while (_peek(rest) != ']')
            {
                int d = _get(rest);
                if (d >= '0' && d <= '9')
                {
                    index = index * 10 + (d - '0');
                }
                else
                {
                    return SLANG_E_INVALID_ARG;
                }
            }

            if (_peek(rest) != ']')
                return SLANG_E_INVALID_ARG;
            _get(rest);

            cursor = cursor.getElement(index);
            state = ALLOW_DOT | ALLOW_SUBSCRIPT;
            continue;
        }
        else
        {
            const char* nameBegin = rest;
            for (;;)
            {
                switch (_peek(rest))
                {
                default:
                    _get(rest);
                    continue;

                case -1:
                case '.':
                case '[':
                    break;
                }
                break;
            }
            char const* nameEnd = rest;
            ShaderCursor newCursor;
            cursor.getField(nameBegin, nameEnd, newCursor);
            cursor = newCursor;
            state = ALLOW_DOT | ALLOW_SUBSCRIPT;
            continue;
        }
    }

    ioCursor = cursor;
    return SLANG_OK;
}

} // namespace fcpw