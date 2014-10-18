#include "REAL.H"
#include "Vector.H"
#include <map>
class LevelGeometry;
template<class T> class LevelData;
class Box;
class FArrayBox;
class FluxBox;
class StructurePool;
class PhysBCUtil;


// Computes N = sqrt(max|db/dz|) over the entire level.
Real computeMaxBVFreq (const LevelData<FluxBox>& a_buoyancy,
                       const LevelGeometry&      a_levGeo,
                       const Real                a_pow);


// Computes sqrt(-db/dz)^p over one grid.
void computeBVFreq (FArrayBox&           a_N,
                    const FluxBox&       a_B,
                    const Box&           a_destBox,
                    const LevelGeometry& a_levGeo,
                    const Real           a_pow);

// Computes sqrt(-db/dz)^p over the entire level.
void computeBVFreq (LevelData<FArrayBox>&     a_N,
                    const LevelData<FluxBox>& a_B,
                    const LevelGeometry&      a_levGeo,
                    const Real                a_pow);

// Finds the max speed, c0 (smallest Nsq/c^2) and the corresponding structure
// function, phi0 (corresponding eigenfunction) of phi'' + (Nsq/c^2) phi = 0.
// This level's grids must span the vertical domain.
void solveVertEigenProblem (LevelData<FArrayBox>&       a_c0,
                            LevelData<FArrayBox>&       a_phi0,
                            const LevelData<FArrayBox>& a_Nsq,
                            const LevelGeometry&        a_levGeo);



// These functions were just copied from old code....

class StructurePool;

// Computes N = sqrt(-db/dz) over one grid.
// This is a more general version...N and B can have (almost) any centering.
// TODO: Make a general mapped version that takes the place of all BV funcs.
void computeBVFreqCartesian (FArrayBox&           a_N,
                             const FArrayBox&     a_B,
                             const Box&           a_destBox,
                             const LevelGeometry& a_levGeo,
                             const Real           a_pow);

// Finds the max speed, c0 (smallest eigenvalue) and the corresponding structure
// function, phi0 (corresponding eigenfunction) of phi'' + (Nsq/c^2) phi = 0.
void solveVertEigenProblemCartesian (Real&            a_c0,
                                     FArrayBox&       a_phi0,
                                     const FArrayBox& a_Nsq,
                                     const Box&       a_valid,
                                     const Real       a_dz);

// Fills a_A0 with the horizontal sech^2 structure.
void fillHorizontalStructure (FArrayBox&           a_A0,
                              StructurePool&       a_pool,
                              const LevelGeometry& a_levGeo,
                              const Real           a_xcenter,
                              const Real           a_amp0);


// -----------------------------------------------------------------------------
// This class retrieves the results of the vertical structure eigenproblems
// for each level and ensures that the equations are only solved one time.
// This is important because the eigenproblem is not easily scalable.
// -----------------------------------------------------------------------------
class StructurePool
{
    // This is what will be stored in the pool - the eigenproblem results.
    typedef struct {
        Real c0, alpha, beta, gamma;
        FArrayBox* CCNsqPtr;
        FArrayBox* CCphiPtr;
    } ElementType;

public:
    // Constructor
    StructurePool ();

    // Destructor
    virtual ~StructurePool ();

    // Sets pointers to objects needed by createNewLevel. If there is any chance
    // that createNewLevel will be triggered, be sure to call this function
    // first or you'll be sorry.
    virtual void setGeometry (const LevelGeometry* a_levGeoPtr,
                              const PhysBCUtil*    a_physBCPtr);

    // Basic accessors
    virtual const LevelGeometry& getLevGeo () const       {return *m_levGeoPtr;}
    virtual const PhysBCUtil& getPhysBCUtil () const      {return *m_physBCPtr;}

    virtual Real c0 (const Box& a_box)                    {return this->getStructureData(a_box).c0;}
    virtual Real c0 (const Box& a_box) const              {return this->getStructureData(a_box).c0;}

    virtual Real alpha (const Box& a_box)                 {return this->getStructureData(a_box).alpha;}
    virtual Real alpha (const Box& a_box) const           {return this->getStructureData(a_box).alpha;}

    virtual Real beta (const Box& a_box)                  {return this->getStructureData(a_box).beta;}
    virtual Real beta (const Box& a_box) const            {return this->getStructureData(a_box).beta;}

    virtual Real gamma (const Box& a_box)                 {return this->getStructureData(a_box).gamma;}
    virtual Real gamma (const Box& a_box) const           {return this->getStructureData(a_box).gamma;}

    virtual const FArrayBox& Nsq (const Box& a_box)       {return *(this->getStructureData(a_box).CCNsqPtr);}
    virtual const FArrayBox& Nsq (const Box& a_box) const {return *(this->getStructureData(a_box).CCNsqPtr);}

    virtual const FArrayBox& phi (const Box& a_box)       {return *(this->getStructureData(a_box).CCphiPtr);}
    virtual const FArrayBox& phi (const Box& a_box) const {return *(this->getStructureData(a_box).CCphiPtr);}

protected:
    // The data accessor for internal use. Checks if a call to createNewLevel
    // is needed and performs some sanity checks on the pool in debug mode.
    virtual ElementType getStructureData (const Box& a_box);
    virtual const ElementType getStructureData (const Box& a_box) const;

    // This is called if a level's data is not yet in the pool. If there is any
    // chance that this function will be triggered, be sure to call setGeometry
    // first or you'll be sorry. This is the workhorse of the class - it solves
    // the eigenprobem using LAPACK.
    virtual ElementType createNewLevel (const Box& a_box);


    // Member variables --------------------------------------------------------

    // The data pool. This links a level's domainBox with its data.
    std::map<Box, ElementType> m_dataMap;

    // A record of all pointers that we will need to delete.
    Vector<ElementType>        m_deletePtrs;

    // These are needed by createNewLevel.
    const LevelGeometry*       m_levGeoPtr;
    const PhysBCUtil*          m_physBCPtr;
};
