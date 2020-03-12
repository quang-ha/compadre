#ifndef _COMPADRE_GMLS_POISSONPURENEUMANN_SOURCES_HPP_
#define _COMPADRE_GMLS_POISSONPURENEUMANN_SOURCES_HPP_

#include <Compadre_SourcesT.hpp>
#include <Compadre_XyzVector.hpp>

namespace Compadre {

class GMLS_PoissonPureNeumannPhysics;

class GMLS_PoissonPureNeumannSources : public SourcesT{
    protected:
        typedef Compadre::ParticlesT particle_type;
        GMLS_PoissonPureNeumannPhysics* _physics;

    public:
        GMLS_PoissonPureNeumannSources(Teuchos::RCP<particle_type> particles, mvec_type* b = NULL) :
                SourcesT(particles, b) {};

        virtual ~GMLS_PoissonPureNeumannSources() {};

        virtual void evaluateRHS(local_index_type field_one, local_index_type field_two = -1, scalar_type time = 0.0);

        virtual std::vector<InteractingFields> gatherFieldInteractions();

        void setPhysics(Teuchos::RCP<GMLS_PoissonPureNeumannPhysics> physics) { _physics = physics.getRawPtr();}
};

}

#endif
