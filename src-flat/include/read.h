
#ifndef __READ_H
#define __READ_H

simflat_t *fromFileASCII(command_t cmd, struct pmd_base_potential_t *pot);
simflat_t *fromFileGzip(command_t cmd, struct pmd_base_potential_t *pot);
simflat_t *fromFileTim(command_t cmd, struct pmd_base_potential_t *pot);

#endif
