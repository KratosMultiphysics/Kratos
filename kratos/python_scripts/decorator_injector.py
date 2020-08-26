import inspect

def InjectIntoAllContainer(container, function):
    InjectIntoContainer(container, function, lambda *x: True)

def InjectIntoContainer(container, function, expression):
    for symbol_name in container.__dict__:
        symbol = getattr(container, symbol_name)

        print(symbol_name, symbol, callable(symbol), inspect.isclass(symbol))
        if callable(symbol) and not inspect.isclass(symbol):
            if expression(symbol_name):
                if hasattr(symbol, "__wrapped__"):
                    print("Warning: symbol", symbol_name, "is already decorated. Injector may fail.")
                
                setattr(container, symbol_name, function(symbol))

        if callable(symbol) and inspect.isclass(symbol):
            InjectIntoContainer(symbol, function, expression)