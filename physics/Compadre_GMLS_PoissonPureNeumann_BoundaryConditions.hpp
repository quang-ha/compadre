#ifndef _COMPADRE_GMLS_POISSONPURENEUMANN_BOUNDARYCONDITIONS_
#define _COMPADRE_GMLS_POISSONPURENEUMANN_BOUNDARYCONDITIONS_

#include <Compadre_BoundaryConditionsT.hpp>

namespace Compadre {

class ParticlesT;

class GMLS_PoissonPureNeumannBoundaryConditions : public BoundaryConditionsT {
    protected:
        typedef Compadre::ParticlesT particle_type;

    public:
        GMLS_PoissonPureNeumannBoundaryConditions(Teuchos::RCP<particle_type> particles, mvec_type* b = NULL) :
                BoundaryConditionsT(particles, b) {}

        virtual ~GMLS_PoissonPureNeumannBoundaryConditions() {}

        virtual void flagBoundaries();

        virtual void applyBoundaries(local_index_type field_one, local_index_type field_two=-1, scalar_type time=0.0);

        virtual std::vector<InteractingFields> gatherFieldInteractions();
};

}

#endif
