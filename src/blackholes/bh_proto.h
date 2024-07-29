/*black hole functions*/

void reallocate_memory_maxpartbh(void);
void domain_resize_storage_bh(int count_get_bh);

void bh_in(void);
void bh_out(void);

void kernel(double u, double z, double hinv3, double hinv4, double *wk, double *dwk);
void bh_density(void);
void bh_ngb_feedback(void);
void bh_feedback(void);
