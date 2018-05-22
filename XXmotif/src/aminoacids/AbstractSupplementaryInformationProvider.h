/**
 * Abstract base class for all methods to include supplementary
 * information.
 *
 * Each concrete class has to implement the two methods
 * for P-value calculation, 1) for the kmer phase, 2) for the PWM phase.
 */

#ifndef _ABSTRACT_SUPPLEMENTARY_INFORMATION_PROVIDER_H
#define _ABSTRACT_SUPPLEMENTARY_INFORMATION_PROVIDER_H

#include "../refinementPhase/Motif.h"
class AbstractKmer;

class AbstractSupplementaryInformationProvider {

public:

  /* Calculate NN prediction P-value for a motif in kmer representation during
   * the initial extension phase.
   *
   */
	virtual double getPvalue(AbstractKmer const* const kmer, const int &seq,
			const int &startPos) const = 0;

  /* Calculate NN prediction P-value for a motif in PWM representation during
   * the final optimization phase.
   */
	virtual double getPvalue(motif_columns_type const* const cols, const int &seq,
			const int &startPos) const = 0;

};

#endif
