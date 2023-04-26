# Basic Data Tytpes in C++

**basic(atomic) data types** of a programming langauge are the types directly provided by the language, as opposed to the data types defined by the programmer using the mechanisms of the langauge.

Some high-level programming langauges provide basic numerical data types that correspond to the ideal mathematical types $\mathbb{R}$ and $\mathbb{Z}$. But C++ provides only low level basic data types that are directly represented by the computer and directly operated on by the ALU. This data types are `int`, `float` and `double`  correspong to $\tilde{\mathbb{Z}}$ and $\mathbb{F}$, respectively. 


| Type           | Range                    | Implements       | Represents         |
| :------------- | :----------------------- | :--------------- | :----------------- |
| `int`          | \[-2^31^, 2^31^-1\]      | IEEE int         | $\mathbb{Z}$       |
| `unsigned int` | \[0, 2^32^ - 1\]         | -                | $\mathbb{N}$       |
| `float`        | \[-3.4e38, 3.4e38]       | IEEE float       | $\mathbb{R}$       |
| `double`       | \[-1.80e+308, 1.80e+308] | IEEE double      | $\mathbb{R}$       |
| `char`         | ASCII characters         | ASCII characters | letters and others |
| `string`       | strings of ASCII         | -                | -                  |

Table: List of Basic Data Types in C++

```{index} basic data type
```