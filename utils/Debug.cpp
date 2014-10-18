#include <iomanip>
#include "Debug.H"
#include "Printing.H"
#include "FluxBox.H"

#ifndef NDEBUG
#   include "Constants.H"
#   include "SetValLevel.H"
#endif

#ifdef CH_MULTIDIM
	using Chombo::pout;
#endif


// -----------------------------------------------------------------------------
// Puts a checkpoint into pout.*
// -----------------------------------------------------------------------------
void _WriteCheckpoint(std::string   a_func,
                      std::string   a_file,
                      unsigned int  a_line,
                      std::string   a_mark) {
#   ifndef NDEBUG
	    using namespace std;
	    string::size_type	start, end, center;

	    center = a_func.find_first_of("::");
	    start = a_func.rfind(" ", center-1) + 1;
	    end   = a_func.find_first_of("(", center+1) - 1;
	    a_func = a_func.substr(start, end - start + 1) + "()";

	    start = a_file.find_last_of("/") + 1;
	    if(start != string::npos && start < a_file.length())
		    a_file = a_file.substr(start, a_file.length() - start);

	    if(a_mark.length() > 0) pout() << setfill('.');
	    pout() << "CHECK:\t" << setiosflags(ios::left) << setw(80) << a_func << "\t" << setw(30) << a_file;
	    if(a_mark.length() > 0) pout() << setfill(' ');
	    pout() << "[" << setw(4) << a_line << "] " << a_mark << endl;
#   endif // End debug code
}


#ifndef NDEBUG
    // -----------------------------------------------------------------------------
    // Used to make TODO() output readable.
    // -----------------------------------------------------------------------------
    static std::string stripFuncName(const std::string& a_func) {
        std::string::size_type	start, end, center;
        center = a_func.find_last_of("::");
        start = a_func.rfind(" ", center-1) + 1;
        end   = a_func.find_first_of("(", center+1) - 1;
        return a_func.substr(start, end - start + 1) + "()";
    }


    // -----------------------------------------------------------------------------
    // Used to make TODO() output readable.
    // -----------------------------------------------------------------------------
    static std::string stripFileName(const std::string& a_file) {
        std::string::size_type start = a_file.find_last_of("/") + 1;
        if(start != std::string::npos && start < a_file.length())
	        return a_file.substr(start, a_file.length() - start);
        return a_file;
    }
#endif // End debug code



// -----------------------------------------------------------------------------
// Send a todo message to stdout
// -----------------------------------------------------------------------------
void _todo(const char* a_filename, const char* a_funcname, const int a_linenumber) {
#   ifndef NDEBUG
	    static bool wasWarned = false;
	    if(!wasWarned && procID() == 0) {
		    std::cout << color::brown << "TODO:"
		              << color::none << " Implement "
		              << stripFuncName(a_funcname)
		              << " in " << stripFileName(a_filename)
		              << " [" << a_linenumber << "]" << std::endl;
//		    wasWarned = true;
	    }
#   endif // End debug code
}



// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox.
// -----------------------------------------------------------------------------
void checkForNAN (const FArrayBox& a_data, const Box& a_testBox, naninfo_type* a_info)
{
#   ifndef NDEBUG
        CH_assert(a_data.box().sameType(a_testBox));

        if (a_info != NULL) {
            a_info->problemFound = false;
        }
        const int proc = procID();

        Box region = a_testBox;
        region &= a_data.box();

        for (int c = 0; c < a_data.nComp(); ++c) {
            BoxIterator bit(region);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& i = bit();
                Real v = a_data(i,c);

                if (v != v || v <= -1.0e10  || 1.0e10 <= v) {
                    ostringstream str;
                    str << "\n\n" << v <<" found inside " << region << " on proc " << proc << " at comp = " << c << ", index = " << i;
                    pout() << str.str().c_str() << std::endl;

                    if (a_info == NULL) {
                        MayDay::Error(str.str().c_str());
                    } else {
                        a_info->problemFound = true;
                        a_info->pos = i;
                        a_info->comp = c;
                        a_info->box = region;
                        a_info->val = v;
                        a_info->msg = str.str();
                        return;
                    }
                }
            }
        }
#   endif
}


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox.
// -----------------------------------------------------------------------------
void checkForNAN (const FluxBox& a_data, const Box& a_testBox, naninfo_type* a_info)
{
#   ifndef NDEBUG
        CH_assert(a_testBox.type() == IntVect::Zero);

        if (a_info != NULL) {
            a_info->problemFound = false;
        }
        const int proc = procID();

        Box CCregion = a_testBox;
        CCregion &= a_data.box();

        for (int d = 0; d < SpaceDim; ++d) {
            Box region = CCregion;
            region.surroundingNodes(d);
            region &= a_data[d].box();

            for (int c = 0; c < a_data.nComp(); ++c) {
                BoxIterator bit(region);
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& i = bit();
                    Real v = a_data[d](i,c);

                    if (v != v || v <= -1.0e10  || 1.0e10 <= v) {
                        ostringstream str;
                        str << "\n\n" << v <<" found inside " << region << " on proc " << proc << " at comp = " << c << ", index = " << i;
                        pout() << str.str().c_str() << std::endl;

                        if (a_info == NULL) {
                            MayDay::Error(str.str().c_str());
                        } else {
                            a_info->problemFound = true;
                            a_info->pos = i;
                            a_info->comp = c;
                            a_info->box = region;
                            a_info->val = v;
                            a_info->msg = str.str();
                            return;
                        }
                    }
                }
            }
        }
#   endif
}

// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the valid regions.
// -----------------------------------------------------------------------------
void checkForValidNAN (const LevelData<FArrayBox>& a_data)
{
#   ifndef NDEBUG
        const DisjointBoxLayout& grids = a_data.getBoxes();
        DataIterator dit = grids.dataIterator();

        naninfo_type info;
        info.problemFound = false;

        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& d = a_data[dit];
            Box valid = grids[dit];
            valid &= d.box();

            checkForNAN(d, valid, &info);
            if (info.problemFound == true) break;
        }

        Vector<int> vpf(numProc(), 0);
        int pf = info.problemFound? 1: 0;
        const int srcProc = uniqueProc(SerialTask::compute);
        gather(vpf, pf, srcProc);
        if (procID() == srcProc) {
            for (int idx = 0; idx < vpf.size(); ++idx) {
                if (vpf[idx] == 1) {
                    pf = 1;
                    break;
                }
            }
        }
        broadcast(pf, srcProc);

        if (pf == 1) {
            const LevelData<FArrayBox>& NANData = a_data;
            writeLevelHDF5(NANData, 0.0, false);

            // const LevelData<FArrayBox>& NANDataOneGhost = a_data;
            // writeLevelHDF5(NANDataOneGhost, 0.0, true);

            MayDay::Error(info.msg.c_str());
        }
#   endif
}


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the valid regions.
// -----------------------------------------------------------------------------
void checkForValidNAN (const LevelData<FluxBox>& a_data)
{
#   ifndef NDEBUG
        const DisjointBoxLayout& grids = a_data.getBoxes();
        DataIterator dit = grids.dataIterator();

        naninfo_type info;
        info.problemFound = false;

        for (dit.reset(); dit.ok(); ++dit) {
            const FluxBox& d = a_data[dit];
            Box valid = grids[dit];
            valid &= d.box();

            checkForNAN(d, valid, &info);
            if (info.problemFound == true) break;
        }

        Vector<int> vpf(numProc(), 0);
        int pf = info.problemFound? 1: 0;
        const int srcProc = uniqueProc(SerialTask::compute);
        gather(vpf, pf, srcProc);
        if (procID() == srcProc) {
            for (int idx = 0; idx < vpf.size(); ++idx) {
                if (vpf[idx] == 1) {
                    pf = 1;
                    break;
                }
            }
        }
        broadcast(pf, srcProc);

        if (pf == 1) {
            const LevelData<FluxBox>& NANData = a_data;
            writeLevelHDF5(NANData, 0.0, false);

            // const LevelData<FluxBox>& NANDataOneGhost = a_data;
            // writeLevelHDF5(NANDataOneGhost, 0.0, true);

            MayDay::Error(info.msg.c_str());
        }
#   endif
}


// -----------------------------------------------------------------------------
// Initializes data holders to NAN
// -----------------------------------------------------------------------------
#ifndef NDEBUG
    // Set a_levels[a_min:a_max] to NAN.
    void debugInitLevels (Vector<BoxLayoutData<FArrayBox>* >& a_levels,
                          int                                 a_min,
                          int                                 a_max)
    {
        setValLevels(a_levels, a_min, a_max, quietNAN);
    }

    // Set a_level to NAN.
    void debugInitLevel (BoxLayoutData<FArrayBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }

    // Set a_level to NAN.
    void debugInitLevel (BoxLayoutData<FluxBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }

    // Set FArrayBox to NAN.
    void debugInit (FArrayBox& a_fab)
    {
        a_fab.setVal(quietNAN);
    }

    // Set FluxBox to NAN.
    void debugInit (FluxBox& a_flub)
    {
        a_flub.setVal(quietNAN);
    }
#endif
