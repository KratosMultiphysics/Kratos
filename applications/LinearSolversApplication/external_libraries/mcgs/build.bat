if [%1]==[] (
    set buildType=Release
) else (
    set buildType=%1
)

cmake -H. -Bbuild -DCMAKE_BUILD_TYPE:STRING=%buildType% -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON -DMCGS_BUILD_TESTS:BOOL=ON
cmake --build build --config %buildType%
