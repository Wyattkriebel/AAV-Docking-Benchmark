#include <cstdio>
#include <cstdlib>
#include "molfile_plugin.h"
extern "C" int molfile_pdbxplugin_init(void);
extern "C" int molfile_pdbxplugin_register(void *, vmdplugin_register_cb);
extern "C" int molfile_pdbxplugin_fini(void);
static molfile_plugin_t *api = NULL;
static int callback(void *, vmdplugin_t *p) { api = (molfile_plugin_t *)p; return 0; }
int main(int argc, char **argv) {
    if (argc != 2) return 2;
    molfile_pdbxplugin_init(); molfile_pdbxplugin_register(NULL, callback);
    int n = 0, flags = 0;
    void *file = api->open_file_read(argv[1], "cif", &n);
    if (file == NULL || n <= 0) return 3;
    molfile_atom_t *atoms = (molfile_atom_t *)calloc(n, sizeof(molfile_atom_t));
    molfile_timestep_t timestep = {};
    timestep.coords = (float *)calloc(3*n, sizeof(float));
    int structure = api->read_structure(file, &flags, atoms);
    int coordinates = api->read_next_timestep(file, n, &timestep);
    if (structure != MOLFILE_SUCCESS || coordinates != MOLFILE_SUCCESS) return 4;
    printf("READ_COUNT\t%d\n", n);
    for (int i=0; i<n; ++i)
        printf("ATOM\t%s\t%d\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\n",
            atoms[i].chain_auth, atoms[i].resid_auth, atoms[i].insertion,
            atoms[i].resname, atoms[i].name, atoms[i].atomicnumber,
            timestep.coords[3*i], timestep.coords[3*i+1], timestep.coords[3*i+2]);
    api->close_file_read(file); molfile_pdbxplugin_fini();
    free(atoms); free(timestep.coords); return 0;
}
