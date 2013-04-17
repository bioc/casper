void addRead2Frag(const char *qname, const char *chr, int start, int strand, const char *cigar, int totF, read_t *frags, int read);

int *procCigar(char *cigar, int *cigs);
