def InjectIntoAllModule(module, function):
    InjectIntoModule(module, function, lambda *x: True)

def InjectIntoModule(module, function, expression):
    for symbol_name in module.__dict__:
        symbol = getattr(module, symbol_name)

        if callable(symbol):
            if expression(symbol_name):
                is_wrapped = hasattr(symbol, "__wrapped__")

                if is_wrapped:
                    print("Warning: function", symbol_name, "is laready decorated. Injectior may fail.")
                
                setattr(module, symbol_name, function(symbol))