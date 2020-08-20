def InjectIntoModule(module, function, expression):
    for symbol_name in module.__dict__:
        symbol = getattr(module, symbol_name)
        if callable(symbol):
            if expression(symbol_name):
                setattr(module, symbol_name, function(symbol))