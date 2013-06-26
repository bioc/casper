typedef struct read_t {
  //int len;
  //int st_1;
  //int end_1;
  int len_1;
  int strand_1;
  //char *cigar_1;
  //int cigar_1;
  //char *chr_1;
  //int st_2;
  //int end_2;
  int len_2;
  int strand_2;
  //char *cigar_2;
  //int cigar_2;
  //char *chr_2;
  //char *qname;
  int nreads;
} read_t;

typedef struct path_t {
  int len;
  char *path;
  int *exons;
  int *starts;
  int *rids;
  int *tx;
  int *rnk;
  int pos;
  int nexon;
  int count;
  int incomp;
  int island;
} path_t;

int verbose;
typedef struct var_t {
  int nex;
  double exp;
  int len;
  int strand;
  int *exst;
  int *exen;
  int *exid;
} var_t;

typedef struct gene_t {
  int nvar;
  double exp;
  var_t *vars;
  char *chr;
} gene_t;
