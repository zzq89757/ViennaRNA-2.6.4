#ifndef VIENNA_RNA_PACKAGE_PAIR_MAT_H
#define VIENNA_RNA_PACKAGE_PAIR_MAT_H

#include <ctype.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/fold_vars.h>

#define NBASES 8
/*@notnull@*/

#ifndef INLINE
# ifdef __GNUC__
#  define INLINE inline
# else
#  define INLINE
# endif
#endif

static const char Law_and_Order[]         = "_ACGUTXKI";
static int        BP_pair[NBASES][NBASES] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

#define MAXALPHA 20       /* maximal length of alphabet */

static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
static int    rtype[8] = {
  0, 2, 1, 4, 3, 6, 5, 7
};

#ifdef _OPENMP
#pragma omp threadprivate(Law_and_Order, BP_pair, alias, pair, rtype)
#endif

/* for backward compatibility */
#define ENCODE(c) encode_char(c)

/*
将碱基字符（例如 A, C, G, U/T 等）转换为对应的数字编码
*/
static INLINE int
encode_char(char c)
{
  /* return numerical representation of base used e.g. in pair[][] */
  int code;

  c = toupper(c);
  /*
  如果 energy_set 大于 0，使用简单的字符编码方案。
  将字符变量c 转换为相对于字母 'A' 的偏移量，并加 1。
  例如，'A' -> 1, 'B' -> 2, ..., 'Z' -> 26。
  */
  if (energy_set > 0) {
    code = (int)(c - 'A') + 1;
  } else {
    /*
    如果 energy_set 小于或等于 0，使用自定义的编码方案。
    Law_and_Order[10] = "_ACGUTXKI"
    如果字符串变量c 未找到，编码为 0；否则计算字符在 Law_and_Order 中的位置。
    */
    const char *pos;
    // strchr(Law_and_Order, c) 查找字符 c 在 Law_and_Order 中的位置。
    // 如果找到，返回指向该位置的指针；否则返回 NULL。
    pos = strchr(Law_and_Order, c);
    if (pos == NULL)
      code = 0;
    else
      // pos - Law_and_Order:计算指针之间的距离，即 pos 相对于 Law_and_Order 起始位置的偏移量（索引）
      code = (int)(pos - Law_and_Order);
    // 如果编码大于 5，设为 0。
    if (code > 5)
      code = 0;
    // 如果编码大于 4，减少 1，以使 'T' 和 'U' 等效（即将 'T' 和 'U' 都编码为 4）。
    if (code > 4)
      code--;           /* make T and U equivalent */
  }

  return code;
}


/*@+boolint +charint@*/
/*@null@*/
extern char *nonstandards;


/*
  配对矩阵
*/
static INLINE void
make_pair_matrix(void)
{

  /* i, j：循环变量。
  alias[]：碱基的别名，用于处理非标准碱基。
  pair[][]：碱基对配对矩阵。
  rtype[]：反向配对矩阵。
  energy_set：指定的能量集。
  nonstandards：允许的非标准碱基对。
  NBASES：碱基的总数:8
  MAXALPHA：最大字母数量:20（根据能量集可能有所不同）。
  BP_pair[][]：默认的碱基对配对矩阵。
  */

  int i, j;
  // 默认能量集 (energy_set == 0)
  if (energy_set == 0) {
    // 将前5个碱基设置为其本身。
    // 将非标准碱基（例如X和K）设置为标准碱基的别名。
    for (i = 0; i < 5; i++)
      alias[i] = (short)i;
    alias[5]  = 3;  /* X <-> G */
    alias[6]  = 2;  /* K <-> C */
    alias[7]  = 0;  /* I <-> default base '@' */
    for (i = 0; i < NBASES; i++)
      for (j = 0; j < NBASES; j++)
        pair[i][j] = BP_pair[i][j];
    // 如果noGU为真，禁止G-U配对。
    if (noGU)
      pair[3][4] = pair[4][3] = 0;
    // 如果nonstandards不为空，允许非标准碱基对
    if (nonstandards != NULL) {
      /* allow nonstandard bp's */
      for (i = 0; i < (int)strlen(nonstandards); i += 2)
        pair[encode_char(nonstandards[i])]
        [encode_char(nonstandards[i + 1])] = 7;
    }
    // 设置反向配对矩阵
    for (i = 0; i < NBASES; i++)
      for (j = 0; j < NBASES; j++)
        rtype[pair[i][j]] = pair[j][i];
  } else {
    for (i = 0; i <= MAXALPHA; i++)
      for (j = 0; j <= MAXALPHA; j++)
        pair[i][j] = 0;
    // 非默认能量集
    if (energy_set == 1) {
      // 设置碱基别名：A <-> G，B <-> C。
      // 初始化配对矩阵：AB <-> GC，BA <-> CG
      for (i = 1; i < MAXALPHA; ) {
        alias[i++]  = 3;      /* A <-> G */
        alias[i++]  = 2;      /* B <-> C */
      }
      for (i = 1; i < MAXALPHA; i++) {
        pair[i][i + 1] = 2;       /* AB <-> GC */
        i++;
        pair[i][i - 1] = 1;       /* BA <-> CG */
      }
    } else if (energy_set == 2) {
      // 设置碱基别名：A <-> A，B <-> U。
      // 初始化配对矩阵：AB <-> AU，BA <-> UA。
      for (i = 1; i < MAXALPHA; ) {
        alias[i++]  = 1;      /* A <-> A*/
        alias[i++]  = 4;      /* B <-> U */
      }
      for (i = 1; i < MAXALPHA; i++) {
        pair[i][i + 1] = 5;       /* AB <-> AU */
        i++;
        pair[i][i - 1] = 6;       /* BA <-> UA */
      }
    } else if (energy_set == 3) {
      // 设置碱基别名：A <-> G，B <-> C，C <-> A，D <-> U。
      // 初始化配对矩阵：AB <-> GC，BA <-> CG，CD <-> AU，DC <-> UA
      for (i = 1; i < MAXALPHA - 2; ) {
        alias[i++]  = 3;    /* A <-> G */
        alias[i++]  = 2;    /* B <-> C */
        alias[i++]  = 1;    /* C <-> A */
        alias[i++]  = 4;    /* D <-> U */
      }
      for (i = 1; i < MAXALPHA - 2; i++) {
        pair[i][i + 1] = 2;     /* AB <-> GC */
        i++;
        pair[i][i - 1] = 1;     /* BA <-> CG */
        i++;
        pair[i][i + 1] = 5;     /* CD <-> AU */
        i++;
        pair[i][i - 1] = 6;     /* DC <-> UA */
      }
    } else {
      vrna_message_error("What energy_set are YOU using??");
    }
    //  设置反向配对矩阵 无论是哪种能量集，最后都会设置rtype[pair[i][j]] = pair[j][i]。
    for (i = 0; i <= MAXALPHA; i++)
      for (j = 0; j <= MAXALPHA; j++)
        rtype[pair[i][j]] = pair[j][i];
  }
}

/*
  用 encode_char 函数将每个字符编码为数字，并将结果存储在数组 S 中,how为0时S[0]为序列长度，否则为序列末尾碱基对应的数值
*/
static INLINE short *
encode_sequence(const char  *sequence,
                short       how)
{ 
  /*
  i 和 l：用于循环和序列长度。
  S：分配内存以存储编码后的序列，长度为 l + 2。这里额外添加两个元素的目的是为了处理循环结构中的边界情况。
  */
  unsigned int  i, l = (unsigned int)strlen(sequence);
  short         *S = (short *)vrna_alloc(sizeof(short) * (l + 2));

  switch (how) {
    /* standard encoding as always used for S */
    // 遍历输入序列 sequence，调用 encode_char 函数将每个字符编码为数字，并将结果存储在数组 S 中。
    // S[l + 1] = S[1]：将序列首尾相连，用于处理循环结构。
    // S[0] = (short)l：存储序列的长度在 S[0] 中。
    case 0:
      for (i = 1; i <= l; i++)    /* make numerical encoding of sequence */
        S[i] = (short)encode_char(sequence[i - 1]);
      S[l + 1]  = S[1];
      S[0]      = (short)l;
      break;
    /* encoding for mismatches of nostandard bases (normally used for S1) */
    // 同样遍历输入序列 sequence，但是通过 alias 数组将每个字符的编码映射为新的编码，并将结果存储在数组 S 中。
    // S[l + 1] = S[1]：同样将序列首尾相连。
    // S[0] = S[l]：存储序列末尾元素的值作为 S[0]，用于处理特定情况。
    case 1:
      for (i = 1; i <= l; i++)
        S[i] = alias[(short)encode_char(sequence[i - 1])];
      S[l + 1]  = S[1];
      S[0]      = S[l];
      break;
  }

  return S;
}


#endif /* VIENNA_RNA_PACKAGE_PAIR_MAT_H */
