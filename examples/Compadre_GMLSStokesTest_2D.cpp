#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>

#include <Compadre_GlobalConstants.hpp>
#include <Kokkos_Core.hpp>
#include <Compadre_ProblemT.hpp>
#include <Compadre_ParticlesT.hpp>
#include <Compadre_FieldManager.hpp>
#include <Compadre_EuclideanCoordsT.hpp>
#include <Compadre_NeighborhoodT.hpp>
#include <Compadre_FieldT.hpp>
#include <Compadre_XyzVector.hpp>
#include <Compadre_AnalyticFunctions.hpp>
#include <Compadre_FileIO.hpp>
#include <Compadre_ParameterManager.hpp>

#include <Compadre_GMLS.hpp>

#include <Compadre_GMLS_Stokes_2D_Operator.hpp>
#include <Compadre_GMLS_Stokes_2D_BoundaryConditions.hpp>
#include <Compadre_GMLS_Stokes_2D_Sources.hpp>

typedef int LO;
typedef long GO;
typedef double ST;
typedef Compadre::XyzVector xyz_type;
typedef Compadre::EuclideanCoordsT CT;

int main(int argc, char* args[]) {
#ifdef TRILINOS_LINEAR_SOLVES
    Teuchos::RCP<Compadre::ParameterManager> parameter_manager;
    if (argc > 1) {
        parameter_manager = Teuchos::rcp(new Compadre::ParameterManager(argc, args));
    } else {
        parameter_manager = Teuchos::rcp(new Compadre::ParameterManager());
        std::cout << "WARNING: No parameter list given. Default parameters used." << std::endl;
    }
    if (parameter_manager->helpRequested()) return 0;
    if (parameter_manager->parseError()) return -1;

    Teuchos::RCP<Teuchos::ParameterList> parameters = parameter_manager->getList();

    Teuchos::GlobalMPISession mpi(&argc, &args);
    Teuchos::oblackholestream bstream;
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

    Kokkos::initialize(argc, args);

    Teuchos::RCP<Teuchos::Time> FirstReadTime = Teuchos::TimeMonitor::getNewCounter("1st Read Time");
    Teuchos::RCP<Teuchos::Time> AssemblyTime = Teuchos::TimeMonitor::getNewCounter("Assembly Time");
    Teuchos::RCP<Teuchos::Time> SolvingTime = Teuchos::TimeMonitor::getNewCounter("Solving Time");
    Teuchos::RCP<Teuchos::Time> WriteTime = Teuchos::TimeMonitor::getNewCounter("Write Time");
    Teuchos::RCP<Teuchos::Time> SecondReadTime = Teuchos::TimeMonitor::getNewCounter("2nd Read Time");

    local_index_type input_dim = parameters->get<Teuchos::ParameterList>("io").get<local_index_type>("input dimensions");

    // This proceed setting up the problem so that the parameters will be propagated down into the physics and bcs
    try {
        parameters->get<std::string>("solution type");
    } catch (Teuchos::Exceptions::InvalidParameter & exception) {
        parameters->set<std::string>("solution type", "sine");
    }

    {
        std::vector<std::string> fnames(1);
        std::vector<double> hsize(1);
        std::vector<double> velocity_errors(1);
        std::vector<double> pressure_errors(1);
        const std::string filename_prefix = parameters->get<Teuchos::ParameterList>("io").get<std::string>("input file prefix");
        fnames[0] = filename_prefix + "6.nc";
        hsize[0] = 6;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(parameters->get<int>("loop size")>3, "Only three mesh levels available for this problem.");

        for (LO i=0; i<parameters->get<int>("loop size"); i++) {
            int Porder = parameters ->get<Teuchos::ParameterList>("remap").get<int>("porder");
            double h_size = 2.0/hsize[i];

            std::string testfilename(fnames[i]);

            Teuchos::RCP<Compadre::ParticlesT> particles = Teuchos::rcp(new Compadre::ParticlesT(parameters, comm, input_dim));
            const CT* coords = (CT*)(particles->getCoordsConst());

            // Read in data file
            FirstReadTime->start();
            Compadre::FileManager fm;
            fm.setReader(testfilename, particles);
            fm.read();
            FirstReadTime->stop();

            particles->zoltan2Initialize();
            particles->getFieldManager()->createField(2, "velocity");
            particles->getFieldManager()->createField(1, "pressure");
            particles->getFieldManager()->createField(2, "velocity_exact");
            particles->getFieldManager()->createField(1, "pressure_exact");
            particles->getFieldManager()->createField(1, "lm_pressure");

            ST halo_size;
            {
                if (parameters->get<Teuchos::ParameterList>("halo").get<bool>("dynamic")) {
                    halo_size = h_size*parameters->get<Teuchos::ParameterList>("halo").get<double>("multiplier");
                } else {
                    halo_size = parameters->get<Teuchos::ParameterList>("halo").get<double>("size");
                }
                particles->buildHalo(halo_size);
                particles->createDOFManager();

                particles->createNeighborhood();

                LO neighbors_needed = Compadre::GMLS::getNP(Porder, 2);

                particles->getNeighborhood()->constructAllNeighborLists(particles->getCoordsConst()->getHaloSize(),
                        parameters->get<Teuchos::ParameterList>("neighborhood").get<std::string>("search type"),
                        true /*dry run for sizes*/,
                        neighbors_needed,
                        parameters->get<Teuchos::ParameterList>("neighborhood").get<double>("cutoff multiplier"),
                        parameters->get<Teuchos::ParameterList>("neighborhood").get<double>("size"),
                        parameters->get<Teuchos::ParameterList>("neighborhood").get<bool>("uniform radii"),
                        parameters->get<Teuchos::ParameterList>("neighborhood").get<double>("radii post search scaling"));

                auto max_h = particles->getNeighborhood()->computeMaxHSupportSize(true);
                auto min_neighbors = particles->getNeighborhood()->computeMinNumNeighbors(true);
                if (comm->getRank()==0) {
                    std::cout << "max _h " << max_h << std::endl;
                    std::cout << "min neighbors " << min_neighbors << std::endl;
                }

                // Iterative solver for the problem
                Teuchos::RCP<Compadre::ProblemT> problem = Teuchos::rcp(new Compadre::ProblemT(particles));

                // construct physics, sources and boundary conditions
                Teuchos::RCP<Compadre::GMLS_Stokes2DPhysics> physics =
                    Teuchos::rcp(new Compadre::GMLS_Stokes2DPhysics(particles, Porder));
                Teuchos::RCP<Compadre::GMLS_Stokes2DSources> source =
                    Teuchos::rcp(new Compadre::GMLS_Stokes2DSources(particles));
                Teuchos::RCP<Compadre::GMLS_Stokes2DBoundaryConditions> bcs =
                    Teuchos::rcp(new Compadre::GMLS_Stokes2DBoundaryConditions(particles));

                // set physics, sources and boundary conditions in the problem
                problem->setPhysics(physics);
                source->setPhysics(physics);

                problem->setSources(source);
                problem->setBCS(bcs);

                // assembly
                AssemblyTime->start();
                problem->initialize();
                AssemblyTime->stop();

                // solving
                SolvingTime->start();
                problem->solve();
                SolvingTime->stop();

                // // compute the residual
                // Compadre::CurlCurlPolyTest v_function;
                // Compadre::SecondOrderBasis p_function;
                // particles->getFieldManager()->getFieldByName("velocity")->localInitFromVectorFunction(&v_function);
                // particles->getFieldManager()->getFieldByName("pressure")->localInitFromScalarFunction(&p_function);
                // auto lm_view = particles->getFieldManager()->getFieldByName("lm_pressure")->getMultiVectorPtr()->getLocalView<Compadre::host_view_type>();
                // for (int jj=0; jj<lm_view.extent(0); ++jj) {
                //     lm_view(jj,0) = 0.0;
                // }
                // particles->getFieldManager()->updateFieldsHaloData();
                // problem->residual();
            }

            Teuchos::RCP<Compadre::AnalyticFunction> velocity_function, pressure_function;
            if (parameters->get<std::string>("solution type")=="sine") {
                velocity_function = Teuchos::rcp_static_cast<Compadre::AnalyticFunction>(Teuchos::rcp(new Compadre::CurlCurl2DSineTest));
                pressure_function = Teuchos::rcp_static_cast<Compadre::AnalyticFunction>(Teuchos::rcp(new Compadre::SineProducts(2)));
            } else {
                velocity_function = Teuchos::rcp_static_cast<Compadre::AnalyticFunction>(Teuchos::rcp(new Compadre::CurlCurl2DPolyTest));
                pressure_function = Teuchos::rcp_static_cast<Compadre::AnalyticFunction>(Teuchos::rcp(new Compadre::SecondOrderBasis(2)));
            }

            // In order to comptue the error norm for pure-Neumann pressure field, the mean value of
            // the computed pressure and the exact pressure is required
            ST pressure_val_mean = 0.0;
            ST pressure_exact_mean = 0.0;
            for (int j=0; j<coords->nLocal(); j++) {
                // Obtain the comptued and the exact values
                xyz_type xyz = coords->getLocalCoords(j);
                ST pressure_val = particles->getFieldManagerConst()->getFieldByName("pressure")->getLocalScalarVal(j);
                ST pressure_exact = pressure_function->evalScalar(xyz);
                // Add the up inside each local processor
                pressure_val_mean += pressure_val;
                pressure_exact_mean += pressure_exact;
            }
            // Perform a reduce all to broadcast the sum value to all processors
            ST pressure_global_val_mean, pressure_global_exact_mean;
            Teuchos::Ptr<ST> pressure_global_val_mean_ptr(&pressure_global_val_mean);
            Teuchos::Ptr<ST> pressure_global_exact_mean_ptr(&pressure_global_exact_mean);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, pressure_val_mean, pressure_global_val_mean_ptr);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, pressure_exact_mean, pressure_global_exact_mean_ptr);
            // Take the average on each local processor
            pressure_global_val_mean /= (double)(coords->nGlobalMax());
            pressure_global_exact_mean /= (double)(coords->nGlobalMax());

            // Get some view to write the exact solution to
            Compadre::host_view_type pressure_exact_view = particles->getFieldManager()->getFieldByName("pressure_exact")->getMultiVectorPtr()->getLocalView<Compadre::host_view_type>();
            Compadre::host_view_type velocity_exact_view = particles->getFieldManager()->getFieldByName("velocity_exact")->getMultiVectorPtr()->getLocalView<Compadre::host_view_type>();

            // check solution
            ST error_velocity_norm = 0.0, error_pressure_norm = 0.0;
            // calculate the exact norm
            ST velocity_norm = 0.0, pressure_norm = 0.0;
            // Then calculate the norm for the velocity field and the zero-mean pressure field
            for (int j=0; j<coords->nLocal(); j++) {
                xyz_type xyz = coords->getLocalCoords(j);
                // calculate velocity norm
                std::vector<ST> velocity_computed = particles->getFieldManagerConst()->getFieldByName("velocity")->getLocalVectorVal(j);
                Compadre::XyzVector velocity_exact_xyz = velocity_function->evalVector(xyz);
                std::vector<ST> velocity_exact(2);
                velocity_exact_xyz.convertToStdVector(velocity_exact);
                for (LO id=0; id<2; id++) {
                    error_velocity_norm += (velocity_computed[id] - velocity_exact[id])*(velocity_computed[id] - velocity_exact[id]);
                    velocity_norm += velocity_exact[id]*velocity_exact[id];
                    // Write to output velocity field
                    velocity_exact_view(j, 0) = velocity_exact[0];
                    velocity_exact_view(j, 1) = velocity_exact[1];
                }
                // remove the mean from the pressure field
                ST pressure_val = particles->getFieldManagerConst()->getFieldByName("pressure")->getLocalScalarVal(j) - pressure_global_val_mean;
                ST pressure_exact = pressure_function->evalScalar(xyz) - pressure_global_exact_mean;
                error_pressure_norm += (pressure_exact - pressure_val)*(pressure_exact - pressure_val);
                pressure_norm += pressure_exact*pressure_exact;
                // Write to output pressure field
                pressure_exact_view(j, 0) = pressure_exact;
            }

            // Now sum up and broadcast both the velocity and pressure norm. And also error norm.
            ST velocity_global_norm, pressure_global_norm;
            Teuchos::Ptr<ST> velocity_global_norm_ptr(&velocity_global_norm);
            Teuchos::Ptr<ST> pressure_global_norm_ptr(&pressure_global_norm);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, velocity_norm, velocity_global_norm_ptr);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, pressure_norm, pressure_global_norm_ptr);
            ST error_velocity_global_norm, error_pressure_global_norm;
            Teuchos::Ptr<ST> error_velocity_global_norm_ptr(&error_velocity_global_norm);
            Teuchos::Ptr<ST> error_pressure_global_norm_ptr(&error_pressure_global_norm);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, error_velocity_norm, error_velocity_global_norm_ptr);
            Teuchos::reduceAll<int, ST>(*comm, Teuchos::REDUCE_SUM, error_pressure_norm, error_pressure_global_norm_ptr);

            // Then obtain the global relative error
            error_velocity_global_norm /= velocity_global_norm;
            error_pressure_global_norm /= pressure_global_norm;
            error_velocity_global_norm = sqrt(error_velocity_global_norm);
            error_pressure_global_norm = sqrt(error_pressure_global_norm);
            if (comm->getRank()==0) {
                std::cout << "Global Velocity Norm: " << error_velocity_global_norm << std::endl;
                std::cout << "Global Pressure Norm: " << error_pressure_global_norm << std::endl;
            }
            velocity_errors[i] = error_velocity_global_norm;
            pressure_errors[i] = error_pressure_global_norm;

            WriteTime->start();
            std::string output_filename = parameters->get<Teuchos::ParameterList>("io").get<std::string>("output file prefix") + std::to_string(i) /* loop */ + parameters->get<Teuchos::ParameterList>("io").get<std::string>("output file");
            std::string writetest_output_filename = parameters->get<Teuchos::ParameterList>("io").get<std::string>("output file prefix") + "writetest" + std::to_string(i) /* loop */ + parameters->get<Teuchos::ParameterList>("io").get<std::string>("output file");
            fm.setWriter(output_filename, particles);
            fm.write();
            WriteTime->stop();

            TEUCHOS_TEST_FOR_EXCEPT_MSG(velocity_errors[i]!=velocity_errors[i], "NaN found in error norm.");
            TEUCHOS_TEST_FOR_EXCEPT_MSG(pressure_errors[i]!=pressure_errors[i], "NaN found in error norm.");
            if ((parameters->get<std::string>("solution type")=="sine") || (parameters->get<std::string>("solution type")=="tanh"))  {
                if (i>0) {
                    TEUCHOS_TEST_FOR_EXCEPT_MSG(velocity_errors[i-1]/velocity_errors[i] < 3.5, std::string("Second order not achieved for sine solution of velocity (should be 4). Is: ") + std::to_string(velocity_errors[i-1]/velocity_errors[i]));
                    TEUCHOS_TEST_FOR_EXCEPT_MSG(pressure_errors[i-1]/pressure_errors[i] < 3.5, std::string("Second order not achieved for sine solution of pressure (should be 4). Is: ") + std::to_string(pressure_errors[i-1]/pressure_errors[i]));
                }
            } else {
                TEUCHOS_TEST_FOR_EXCEPT_MSG(velocity_errors[i] > 1e-8, "Second order solution not recovered exactly for velocity.");
                TEUCHOS_TEST_FOR_EXCEPT_MSG(pressure_errors[i] > 1e-8, "Second order solution not recovered exactly for pressure.");
            }
        }
    }
    Teuchos::TimeMonitor::summarize();
    Kokkos::finalize();
#endif
    return 0;
}
