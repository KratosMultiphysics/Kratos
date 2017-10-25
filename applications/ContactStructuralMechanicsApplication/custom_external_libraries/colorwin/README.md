# colorwin
colorwin is a small Windows C++ library under the MIT license that provides scoped coloring of console output.  A single header file with a simple API comprises the library.  The scoped aspect of the library expedites safely adding color to existing console projects.

# colorwin's scoped color control
Scoped coloring means that a colorwin object is created to change the console color and the color is reset back to the original color on destruction of that same object.  Thus, there is less chance of accidentally coloring output you don't mean to.  For example, if your project is using exceptions and using colorwin objects, an exception would reset the console color through C++'s automatic destruction of stack objects (C++ stack unwinding).  If you are adding console coloring to existing projects, you never have to worry about resetting the color back to normal; that means 1 change instead of 2 changes, meaning less opportunities for mistakes and omissions.

# examples

## example 1: DemoAllColors
![example_1_DemoAllColors](https://github.com/jrebacz/colorwin/blob/readme_content/images/example_1_DemoAllColors.png)


```C++
// Show text in each color.
void DemoAllColors()
{
    cout << color(white) << "DemoAllColors()\n";
    cout << color(red) << "red\n";
    cout << color(yellow) << "yellow\n";
    cout << color(green) << "green\n";
    cout << color(cyan) << "cyan\n";
    cout << color(blue) << "blue\n";
    cout << color(magenta) << "magenta\n";
    cout << color(white) << "white\n";
    cout << color(gray) << "gray\n";
    cout << color(grey) << "grey\n";
    cout << color(dark_gray) << "dark_gray\n";
    cout << color(dark_grey) << "dark_grey\n";
}
```

## example 2: DemoException
![example_2_DemoException](https://github.com/jrebacz/colorwin/blob/readme_content/images/example_2_DemoException.png)
```C++
// Demonstrate that scoped colors revert to the original color
// if an exception is thrown "unexpectedly".
void DemoException()
{
    cout << color(white) << "DemoException()\n";
    try
    {
        cout << "This is the original text color.\n";
        withcolor scoped(cyan);
        cout << "Connection " << color(green) << "OK";
        cout << " attempting transfer.\n";
        throw exception("Connection lost"); //oops, I hope the console color doesn't remain cyan!
        cout << "Transfer completed.\n";
    }
    catch (exception &e)
    {
        cout << color(red) << "ERROR: " << e.what() << "\n";
        cout << "Transfer aborted.\n";
    }
    cout << "Text color is back to the original color.\n";
}
```

## example 3: DemoScoped
![example_3_DemoScoped](https://github.com/jrebacz/colorwin/blob/readme_content/images/example_3_DemoScoped.png)
```C++
// Demonstrate multiple levels of scoped console colors.
void DemoScoped()
{
    cout << color(white) << "DemoScoped()\n";
    {
        withcolor scoped(red);
        cout << "|red\n";
        cout << "|red again\n";
        {
            withcolor scoped(yellow);
            cout << "||now yellow\n";
            {
                withcolor scoped(cyan);
                cout << "|||now cyan\n";
                withcolor(white).printf("|||| withcolor(white).printf(...)\n");
                printf("|||::printf cyan\n");

            }
        }
        cout << "|now back to red.\n";
    }
    cout << "now back to normal\n";
}
```

# caveats
* Your console application should protect against multiple threads writing to the console concurrently.  In that case, colorwin should only be used within that protection.
* Dark colors are commented out in the color enumeration to discourage you from using them.  You may comment them in, but they are not very readable.
* Changing the background color is not enabled in colorwin.  Users change the default background color of the command prompt, and Powershell's background is blue by default.  Setting the text background color will result in filled rectangles around text.

