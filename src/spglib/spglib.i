%module spglib
%inline %{

extern void spg_show_symmetry(const double (*lattice)[3], const double (*position)[3], const int* types, const int num_atom, const double symprec);

extern int spg_find_primitive(const double (*lattice)[3], const double (*position)[3], const int* types, const int num_atom, const double symprec);

extern int spg_get_international(char symbol[11], const double lattice[3][3],
                        const double position[][3],
                        const int types[], const int num_atom,
                        const double symprec);

double (*getVectsArray())[3] {
   return (double(*)[3]) malloc(3*3*sizeof(double));
}

void putInDoubleNby3Array(double (*a)[3], int i, int j, double v) {
   a[i][j] = v;
}

double getFromDoubleNby3Array(double (*a)[3], int i, int j) {
   double result = a[i][j];
   return result;
}

void deleteIntArray(int* t) {
   free(t);
}
void deleteDoubleNby3Array(double (*a)[3]) {
   free(a);
}

double (*getPositionsArray(int size))[3] {
   return (double(*)[3]) malloc(size*3*sizeof(double));
}

int* getTypesArray(int size) {
   return (int*) malloc(size*sizeof(int));
}

int getFromIntArray(int* a, int i) {
   return a[i];
}

void putInIntArray(int* a, int i, int v) {
   a[i] = v;
}

%}

