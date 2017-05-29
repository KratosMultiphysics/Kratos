// The MIT License(MIT)
// Copyright(c) 2016  Jeff Rebacz
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files(the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and / or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions :
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef COLORWIN_HPP_INCLUDED
#define COLORWIN_HPP_INCLUDED

#include <Windows.h>
#include <iostream>
#include <stack>

namespace colorwin
{
    // mix colors from wincon.h
    //#define FOREGROUND_BLUE      0x0001 // text color contains blue.
    //#define FOREGROUND_GREEN     0x0002 // text color contains green.
    //#define FOREGROUND_RED       0x0004 // text color contains red.
    //#define FOREGROUND_INTENSITY 0x0008 // text color is intensified.
    enum CW_COLORS
    {
        red = FOREGROUND_RED | FOREGROUND_INTENSITY,
        yellow = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY,
        green = FOREGROUND_GREEN | FOREGROUND_INTENSITY,
        cyan = FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY,
        blue = FOREGROUND_BLUE | FOREGROUND_INTENSITY,
        magenta = FOREGROUND_BLUE | FOREGROUND_RED | FOREGROUND_INTENSITY,
        white = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY,
        gray = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE,
        grey = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE,

        dark_gray = FOREGROUND_INTENSITY,
        dark_grey = FOREGROUND_INTENSITY,

        /* The following dark colors are unreadable, but comment them in and use them if you want. */
        /*
        dark_red = FOREGROUND_RED,
        dark_yellow = FOREGROUND_RED | FOREGROUND_GREEN,
        dark_green = FOREGROUND_GREEN,
        dark_cyan = FOREGROUND_GREEN | FOREGROUND_BLUE,
        dark_blue = FOREGROUND_BLUE,
        dark_magenta = FOREGROUND_BLUE | FOREGROUND_RED,
        */
    };

    // Example usage: std::cout << color(yellow) << "This is a warning color!\n";
    class color
    {
    public:
        color(CW_COLORS color) : m_color(color), m_console_handle(INVALID_HANDLE_VALUE)
        {
            CONSOLE_SCREEN_BUFFER_INFO console_info;
            m_console_handle = GetStdHandle(STD_OUTPUT_HANDLE);
            if (!GetConsoleScreenBufferInfo(m_console_handle, &console_info))
            {
                m_console_handle = GetStdHandle(STD_ERROR_HANDLE);
                if (!GetConsoleScreenBufferInfo(m_console_handle, &console_info)) // maybe standard output device has been redirected, try the standard error device
                {
                    m_console_handle = INVALID_HANDLE_VALUE;
                }
            }
        }

        ~color()
        {
            if (m_console_handle != INVALID_HANDLE_VALUE)
            {
                // Restore the previous color.
                SetConsoleTextAttribute(m_console_handle, get_color_stack().top());
                get_color_stack().pop();
            }
        }

    private:
        void change_color() const
        {
            if (m_console_handle == INVALID_HANDLE_VALUE)
                return; // Can't get console info, can't change color.
            CONSOLE_SCREEN_BUFFER_INFO console_info;
            GetConsoleScreenBufferInfo(m_console_handle, &console_info);
            // save the current attributes for restoration on destruction.
            get_color_stack().push(console_info.wAttributes);
            SetConsoleTextAttribute(m_console_handle, 0x0F & m_color | 0xf0 & console_info.wAttributes);    // save the background color
        }

        color(color &);
        color& operator=(color);

        static std::stack<WORD>& get_color_stack()
        {
            // Use this instead of static member to avoid multiply defined symbols.
            static std::stack<WORD> color_stack;
            return color_stack;
        }

        HANDLE m_console_handle;
        const CW_COLORS m_color;

        friend class withcolor;
        template<typename charT, typename traits> friend std::basic_ostream<charT, traits> & operator<<(std::basic_ostream<charT, traits> &lhs, colorwin::color const &rhs);
    };

    // Example usage : 
    //  {
    //      withcolor scoped(yellow);
    //      cout << "This is a yellow warning!\n";
    //      cout << "This is a second yellow warning!\n";
    //  }
    //  --  or  --
    //      withcolor(yellow).printf("This will be yellow\n");
    class withcolor
    {
    public:
        withcolor(CW_COLORS color) : m_color(color)
        {
            m_color.change_color();
        }

        int printf(const char* format, ...)
        {
            va_list vlist;
            va_start(vlist, format);
            int ret = vprintf(format, vlist);
            va_end(vlist);
            return ret;
        }

        #if _MSC_VER >= 1400    // printf_s offered in Visual Studio 2005
        int printf_s(const char *format, ...)
        {
            va_list vlist;
            va_start(vlist, format);
            int ret = vprintf_s(format, vlist);
            va_end(vlist);
            return ret;
        }
        #endif

    private:
        withcolor(withcolor &);
        withcolor& operator=(withcolor);

        color m_color;
    };

    // cout << color(red) -> operator<<(cout, colorwin::color(red))
    template<typename charT, typename traits> std::basic_ostream<charT, traits> & operator<<(std::basic_ostream<charT, traits> &lhs, colorwin::color const &rhs)
    {
        rhs.change_color();
        return lhs;
    }
}

#endif // COLORWIN_HPP_INCLUDED