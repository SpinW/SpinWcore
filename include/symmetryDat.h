//
// Created by ward_s on 19/04/18.
//

#ifndef SPINW_SYMMETRYDAT_H
#define SPINW_SYMMETRYDAT_H

#include <string>

class ALL_SYM_DAT {
public:
    int symStrLen = 230;
    std::string allStrings[230] = {
            "   1  P 1        : x,y,z",
            "   2  P -1       : -x,-y,-z",
            "   3  P 2        : -x,y,-z",
            "   4  P 21       : -x,y+1/2,-z",
            "   5  C 2        : x+1/2,y+1/2,z; -x,y,-z",
            "   6  P m        : x,-y,z",
            "   7  P c        : x,-y,z+1/2",
            "   8  C m        : x+1/2,y+1/2,z; x,-y,z",
            "   9  C c        : x+1/2,y+1/2,z; x,-y,z+1/2",
            "  10  P 2/m      : -x,y,-z; -x,-y,-z",
            "  11  P 21/m     : -x,y+1/2,-z; -x,-y,-z",
            "  12  C 2/m      : x+1/2,y+1/2,z; -x,y,-z; -x,-y,-z",
            "  13  P 2/c      : -x,y,-z+1/2; -x,-y,-z",
            "  14  P 21/c     : -x,y+1/2,-z+1/2; -x,-y,-z",
            "  15  C 2/c      : x+1/2,y+1/2,z; -x,y,-z+1/2; -x,-y,-z",
            "  16  P 2 2 2    : -x,-y,z; -x,y,-z",
            "  17  P 2 2 21   : -x,-y,z+1/2; -x,y,-z+1/2",
            "  18  P 21 21 2  : -x,-y,z; -x+1/2,y+1/2,-z",
            "  19  P 21 21 21 : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2",
            "  20  C 2 2 21   : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2",
            "  21  C 2 2 2    : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z",
            "  22  F 2 2 2    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z",
            "  23  I 2 2 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z",
            "  24  I 21 21 21 : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2",
            "  25  P m m 2    : -x,-y,z; x,-y,z",
            "  26  P m c 21   : -x,-y,z+1/2; x,-y,z+1/2",
            "  27  P c c 2    : -x,-y,z; x,-y,z+1/2",
            "  28  P m a 2    : -x,-y,z; x+1/2,-y,z",
            "  29  P c a 21   : -x,-y,z+1/2; x+1/2,-y,z",
            "  30  P n c 2    : -x,-y,z; x,-y+1/2,z+1/2",
            "  31  P m n 21   : -x+1/2,-y,z+1/2; x+1/2,-y,z+1/2",
            "  32  P b a 2    : -x,-y,z; x+1/2,-y+1/2,z",
            "  33  P n a 21   : -x,-y,z+1/2; x+1/2,-y+1/2,z",
            "  34  P n n 2    : -x,-y,z; x+1/2,-y+1/2,z+1/2",
            "  35  C m m 2    : x+1/2,y+1/2,z; -x,-y,z; x,-y,z",
            "  36  C m c 21   : x+1/2,y+1/2,z; -x,-y,z+1/2; x,-y,z+1/2",
            "  37  C c c 2    : x+1/2,y+1/2,z; -x,-y,z; x,-y,z+1/2",
            "  38  A m m 2    : x,y+1/2,z+1/2; -x,-y,z; x,-y,z",
            "  39  A b m 2    : x,y+1/2,z+1/2; -x,-y,z; x,-y+1/2,z",
            "  40  A m a 2    : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z",
            "  41  A b a 2    : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z",
            "  42  F m m 2    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x,-y,z",
            "  43  F d d 2    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x+1/4,-y+1/4,z+1/4",
            "  44  I m m 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; x,-y,z",
            "  45  I b a 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z",
            "  46  I m a 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z",
            "  47  P m m m    : -x,-y,z; -x,y,-z; -x,-y,-z",
            "  48  P n n n    : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; -x,-y,-z",
            "  49  P c c m    : -x,-y,z; -x,y,-z+1/2; -x,-y,-z",
            "  50  P b a n    : -x+1/2,-y+1/2,z; -x+1/2,y,-z; -x,-y,-z",
            "  51  P m m a    : -x+1/2,-y,z; -x,y,-z; -x,-y,-z",
            "  52  P n n a    : -x+1/2,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z",
            "  53  P m n a    : -x+1/2,-y,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z",
            "  54  P c c a    : -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z",
            "  55  P b a m    : -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z",
            "  56  P c c n    : -x+1/2,-y+1/2,z; -x,y+1/2,-z+1/2; -x,-y,-z",
            "  57  P b c m    : -x,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z",
            "  58  P n n m    : -x,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z",
            "  59  P m m n    : -x+1/2,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z",
            "  60  P b c n    : -x+1/2,-y+1/2,z+1/2; -x,y,-z+1/2; -x,-y,-z",
            "  61  P b c a    : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z",
            "  62  P n m a    : -x+1/2,-y,z+1/2; -x,y+1/2,-z; -x,-y,-z",
            "  63  C m c m    : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2; -x,-y,-z",
            "  64  C m c a    : x+1/2,y+1/2,z; -x,-y+1/2,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z",
            "  65  C m m m    : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z; -x,-y,-z",
            "  66  C c c m    : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z+1/2; -x,-y,-z",
            "  67  C m m a    : x+1/2,y+1/2,z; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z",
            "  68  C c c a    : x+1/2,y+1/2,z; -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z",
            "  69  F m m m    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z",
            "  70  F d d d    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+3/4,-y+3/4,z; -x+3/4,y,-z+3/4; -x,-y,-z",
            "  71  I m m m    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z",
            "  72  I b a m    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z",
            "  73  I b c a    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z",
            "  74  I m m a    : x+1/2,y+1/2,z+1/2; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z",
            "  75  P 4        : -x,-y,z; -y,x,z",
            "  76  P 41       : -x,-y,z+1/2; -y,x,z+1/4",
            "  77  P 42       : -x,-y,z; -y,x,z+1/2",
            "  78  P 43       : -x,-y,z+1/2; -y,x,z+3/4",
            "  79  I 4        : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z",
            "  80  I 41       : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4",
            "  81  P -4       : -x,-y,z; y,-x,-z",
            "  82  I -4       : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z",
            "  83  P 4/m      : -x,-y,z; -y,x,z; -x,-y,-z",
            "  84  P 42/m     : -x,-y,z; -y,x,z+1/2; -x,-y,-z",
            "  85  P 4/n      : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,-y,-z",
            "  86  P 42/n     : -x+1/2,-y+1/2,z; -y,x+1/2,z+1/2; -x,-y,-z",
            "  87  I 4/m      : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,-y,-z",
            "  88  I 41/a     : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+3/4,x+1/4,z+1/4; -x,-y,-z",
            "  89  P 4 2 2    : -x,-y,z; -y,x,z; -x,y,-z",
            "  90  P 4 21 2   : -x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z",
            "  91  P 41 2 2   : -x,-y,z+1/2; -y,x,z+1/4; -x,y,-z",
            "  92  P 41 21 2  : -x,-y,z+1/2; -y+1/2,x+1/2,z+1/4; -x+1/2,y+1/2,-z+1/4",
            "  93  P 42 2 2   : -x,-y,z; -y,x,z+1/2; -x,y,-z",
            "  94  P 42 21 2  : -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2",
            "  95  P 43 2 2   : -x,-y,z+1/2; -y,x,z+3/4; -x,y,-z",
            "  96  P 43 21 2  : -x,-y,z+1/2; -y+1/2,x+1/2,z+3/4; -x+1/2,y+1/2,-z+3/4",
            "  97  I 4 2 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z",
            "  98  I 41 2 2   : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x+1/2,y,-z+3/4",
            "  99  P 4 m m    : -x,-y,z; -y,x,z; x,-y,z",
            " 100  P 4 b m    : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z",
            " 101  P 42 c m   : -x,-y,z; -y,x,z+1/2; x,-y,z+1/2",
            " 102  P 42 n m   : -x,-y,z; -y+1/2,x+1/2,z+1/2; x+1/2,-y+1/2,z+1/2",
            " 103  P 4 c c    : -x,-y,z; -y,x,z; x,-y,z+1/2",
            " 104  P 4 n c    : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z+1/2",
            " 105  P 42 m c   : -x,-y,z; -y,x,z+1/2; x,-y,z",
            " 106  P 42 b c   : -x,-y,z; -y,x,z+1/2; x+1/2,-y+1/2,z",
            " 107  I 4 m m    : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z ",
            " 108  I 4 c m    : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z+1/2",
            " 109  I 41 m d   : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z",
            " 110  I 41 c d   : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z+1/2",
            " 111  P -4 2 m   : -x,-y,z; y,-x,-z; -x,y,-z",
            " 112  P -4 2 c   : -x,-y,z; y,-x,-z; -x,y,-z+1/2",
            " 113  P -4 21 m  : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z",
            " 114  P -4 21 c  : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z+1/2",
            " 115  P -4 m 2   : -x,-y,z; y,-x,-z; x,-y,z",
            " 116  P -4 c 2   : -x,-y,z; y,-x,-z; x,-y,z+1/2",
            " 117  P -4 b 2   : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z",
            " 118  P -4 n 2   : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z+1/2",
            " 119  I -4 m 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z",
            " 120  I -4 c 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z+1/2",
            " 121  I -4 2 m   : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x,y,-z",
            " 122  I -4 2 d   : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x+1/2,y,-z+3/4",
            " 123  P 4/m m m  : -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z",
            " 124  P 4/m c c  : -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z",
            " 125  P 4/n b m  : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z; -x,-y,-z",
            " 126  P 4/n n c  : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z+1/2; -x,-y,-z",
            " 127  P 4/m b m  : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z; -x,-y,-z",
            " 128  P 4/m n c  : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z",
            " 129  P 4/n m m  : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z; -x,-y,-z",
            " 130  P 4/n c c  : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z+1/2; -x,-y,-z",
            " 131  P 42/m m c : -x,-y,z; -y,x,z+1/2; -x,y,-z; -x,-y,-z",
            " 132  P 42/m c m : -x,-y,z; -y,x,z+1/2; -x,y,-z+1/2; -x,-y,-z",
            " 133  P 42/n b c : -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z; -x,-y,-z",
            " 134  P 42/n n m : -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z",
            " 135  P 42/m b c : -x,-y,z; -y,x,z+1/2; -x+1/2,y+1/2,-z; -x,-y,-z",
            " 136  P 42/m n m : -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2; -x,-y,-z",
            " 137  P 42/n m c : -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z; -x,-y,-z",
            " 138  P 42/n c m : -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z",
            " 139  I 4/m m m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z",
            " 140  I 4/m c m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z",
            " 141  I 41/a m d : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z+1/2; -x,-y,-z",
            " 142  I 41/a c d : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z; -x,-y,-z",
            " 143  P 3        : -y,x-y,z",
            " 144  P 31       : -y,x-y,z+1/3",
            " 145  P 32       : -y,x-y,z+2/3",
            " 146  R 3        : x+1/3,y+2/3,z+2/3; -y,x-y,z",
            " 147  P -3       : -y,x-y,z; -x,-y,-z",
            " 148  R -3       : x+1/3,y+2/3,z+2/3; -y,x-y,z; -x,-y,-z",
            " 149  P 3 1 2    : -y,x-y,z; -y,-x,-z",
            " 150  P 3 2 1    : -y,x-y,z; y,x,-z",
            " 151  P 31 1 2   : -y,x-y,z+1/3; -y,-x,-z+2/3",
            " 152  P 31 2 1   : -y,x-y,z+1/3; y,x,-z",
            " 153  P 32 1 2   : -y,x-y,z+2/3; -y,-x,-z+1/3",
            " 154  P 32 2 1   : -y,x-y,z+2/3; y,x,-z",
            " 155  R 3 2      : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z",
            " 156  P 3 m 1    : -y,x-y,z; -y,-x,z",
            " 157  P 3 1 m    : -y,x-y,z; y,x,z",
            " 158  P 3 c 1    : -y,x-y,z; -y,-x,z+1/2",
            " 159  P 3 1 c    : -y,x-y,z; y,x,z+1/2",
            " 160  R 3 m      : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z",
            " 161  R 3 c      : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z+1/2",
            " 162  P -3 1 m   : -y,x-y,z; -y,-x,-z; -x,-y,-z",
            " 163  P -3 1 c   : -y,x-y,z; -y,-x,-z+1/2; -x,-y,-z",
            " 164  P -3 m 1   : -y,x-y,z; y,x,-z; -x,-y,-z",
            " 165  P -3 c 1   : -y,x-y,z; y,x,-z+1/2; -x,-y,-z",
            " 166  R -3 m     : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z; -x,-y,-z",
            " 167  R -3 c     : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z+1/2; -x,-y,-z",
            " 168  P 6        : -y,x-y,z; -x,-y,z",
            " 169  P 61       : -y,x-y,z+1/3; -x,-y,z+1/2",
            " 170  P 65       : -y,x-y,z+2/3; -x,-y,z+1/2",
            " 171  P 62       : -y,x-y,z+2/3; -x,-y,z",
            " 172  P 64       : -y,x-y,z+1/3; -x,-y,z",
            " 173  P 63       : -y,x-y,z; -x,-y,z+1/2",
            " 174  P -6       : -y,x-y,z; x,y,-z",
            " 175  P 6/m      : -y,x-y,z; -x,-y,z; -x,-y,-z",
            " 176  P 63/m     : -y,x-y,z; -x,-y,z+1/2; -x,-y,-z",
            " 177  P 6 2 2    : -y,x-y,z; -x,-y,z; y,x,-z",
            " 178  P 61 2 2   : -y,x-y,z+1/3; -x,-y,z+1/2; y,x,-z+1/3",
            " 179  P 65 2 2   : -y,x-y,z+2/3; -x,-y,z+1/2; y,x,-z+2/3",
            " 180  P 62 2 2   : -y,x-y,z+2/3; -x,-y,z; y,x,-z+2/3",
            " 181  P 64 2 2   : -y,x-y,z+1/3; -x,-y,z; y,x,-z+1/3",
            " 182  P 63 2 2   : -y,x-y,z; -x,-y,z+1/2; y,x,-z",
            " 183  P 6 m m    : -y,x-y,z; -x,-y,z; -y,-x,z",
            " 184  P 6 c c    : -y,x-y,z; -x,-y,z; -y,-x,z+1/2",
            " 185  P 63 c m   : -y,x-y,z; -x,-y,z+1/2; -y,-x,z+1/2",
            " 186  P 63 m c   : -y,x-y,z; -x,-y,z+1/2; -y,-x,z",
            " 187  P -6 m 2   : -y,x-y,z; x,y,-z; -y,-x,z",
            " 188  P -6 c 2   : -y,x-y,z; x,y,-z+1/2; -y,-x,z+1/2",
            " 189  P -6 2 m   : -y,x-y,z; x,y,-z; y,x,-z ",
            " 190  P -6 2 c   : -y,x-y,z; x,y,-z+1/2; y,x,-z",
            " 191  P 6/m m m  : -y,x-y,z; -x,-y,z; y,x,-z; -x,-y,-z",
            " 192  P 6/m c c  : -y,x-y,z; -x,-y,z; y,x,-z+1/2; -x,-y,-z",
            " 193  P 63/m c m : -y,x-y,z; -x,-y,z+1/2; y,x,-z+1/2; -x,-y,-z",
            " 194  P 63/m m c : -y,x-y,z; -x,-y,z+1/2; y,x,-z; -x,-y,-z",
            " 195  P 2 3      : -x,-y,z; -x,y,-z; z,x,y",
            " 196  F 2 3      : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y",
            " 197  I 2 3      : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y",
            " 198  P 21 3     : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y",
            " 199  I 21 3     : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y",
            " 200  P m -3     : -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z",
            " 201  P n -3     : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; -x,-y,-z",
            " 202  F m -3     : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z",
            " 203  F d -3     : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+1/4,-y+1/4,z; -x+1/4,y,-z+1/4; z,x,y; -x,-y,-z",
            " 204  I m -3     : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z",
            " 205  P a -3     : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z",
            " 206  I a -3     : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z",
            " 207  P 4 3 2    : -x,-y,z; -x,y,-z; z,x,y; y,x,-z",
            " 208  P 42 3 2   : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2",
            " 209  F 4 3 2    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z",
            " 210  F 41 3 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; y+3/4,x+1/4,-z+3/4",
            " 211  I 4 3 2    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z",
            " 212  P 43 3 2   : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+3/4,-z+3/4",
            " 213  P 41 3 2   : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4",
            " 214  I 41 3 2   : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4",
            " 215  P -4 3 m   : -x,-y,z; -x,y,-z; z,x,y; y,x,z",
            " 216  F -4 3 m   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z",
            " 217  I -4 3 m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z",
            " 218  P -4 3 n   : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2",
            " 219  F -4 3 c   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2",
            " 220  I -4 3 d   : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+1/4,z+1/4",
            " 221  P m -3 m   : -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z",
            " 222  P n -3 n   : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y,x,-z+1/2; -x,-y,-z",
            " 223  P m -3 n   : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z",
            " 224  P n -3 m   : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y+1/2,x+1/2,-z; -x,-y,-z",
            " 225  F m -3 m   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z",
            " 226  F m -3 c   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z",
            " 227  F d -3 m   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+3/4,-y+1/4,z+1/2; -x+1/4,y+1/2,-z+3/4; z,x,y; y+3/4,x+1/4,-z+1/2; -x,-y,-z",
            " 228  F d -3 c   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+1/4,-y+3/4,z+1/2; -x+3/4,y+1/2,-z+1/4; z,x,y; y+3/4,x+1/4,-z; -x,-y,-z",
            " 229  I m -3 m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z",
            " 230  I a -3 d   : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4; -x,-y,-z"
    };
};

#endif //SPINW_SYMMETRYDAT_H