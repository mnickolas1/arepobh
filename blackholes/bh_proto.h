/*black hole functions*/

void reallocate_memory_maxpartbh(void);
void domain_resize_storage_bh(int count_get_bh);

void kernel(double u, double hinv3, double hinv4, double *wk, double *dwk);
void bh_density(void);
void bh_ngb_feedback(void);

/*void update_bh_timesteps(void);
void reconstruct_bh_timebins(void);
void update_list_of_active_bh_particles(void);*/
void active_virtual_part_init_alloc(struct ActiveVirtualPart *AVP, const char *name, int *MaxPart);
void active_virtual_part_set(struct ActiveVirtualPart *AVP);
void virtual_part_feedback(void);

void create_particles(void);
void destroy_particles(void);



