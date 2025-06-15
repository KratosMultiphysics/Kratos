import os
import sys
from ui_menu import create_menu


def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(current_dir)
    element_test_path = os.path.join(current_dir, "applications", "GeoMechanicsApplication", "python_scripts", "element_test")
    if element_test_path not in sys.path:
        sys.path.insert(0, element_test_path)

    try:
        create_menu()
    except ImportError as e:
        print(f"[ERROR] Could not run GUI: {e}. Make sure the GeoMechanicsApplication is built and installed correctly.")
        sys.exit(1)

if __name__ == "__main__":
    main()
