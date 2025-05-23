from tkinter import scrolledtext

_log_widget = None

def init_log_widget(widget):
    global _log_widget
    _log_widget = widget

def log_message(msg, level="info"):
    if _log_widget is None:
        print(f"{level.upper()}: {msg}")
        return

    _log_widget.config(state="normal")
    prefix = {"info": "[INFO]", "error": "[ERROR]", "warn": "[WARN]"}
    _log_widget.insert("end", f"{prefix.get(level, '[INFO]')} {msg}\n")
    _log_widget.see("end")
    _log_widget.config(state="disabled")

def clear_log():
    if _log_widget is not None:
        _log_widget.config(state="normal")
        _log_widget.delete("1.0", "end")
        _log_widget.config(state="disabled")

