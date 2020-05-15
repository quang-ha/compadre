#ifndef _COMPADRE_GMLS_STOKES2DSOURCES_HPP_
#define _COMPADRE_GMLS_STOKES2DSOURCES_HPP_

#include <Compadre_SourcesT.hpp>
#include <Compadre_XyzVector.hpp>

namespace Compadre {

class GMLS_Stokes2DPhysics;

class GMLS_Stokes2DSources : public SourcesT{
    protected:
        typedef Compadre::ParticlesT particle_type;
        GMLS_Stokes2DPhysics* _physics;

    public:
        GMLS_Stokes2DSources(Teuchos::RCP<particle_type> particles, mvec_type* b = NULL) :
                SourcesT(particles, b) {};

        virtual ~GMLS_Stokes2DSources() {};

        virtual void evaluateRHS(local_index_type field_one, local_index_type field_two = -1, scalar_type time = 0.0);

        virtual std::vector<InteractingFields> gatherFieldInteractions();

        void setPhysics(Teuchos::RCP<GMLS_Stokes2DPhysics> physics) { _physics = physics.getRawPtr();}
};

}

#endif
