/*
 * FDABase.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <iomanip>
#include <iostream>
#include "FDABase.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "PureInteractionType.h"
#include "Utilities.h"

namespace fda {

template <class Base>
FDABase<Base>::FDABase(ResultType result_type, int syslen, std::string const& result_filename, FDASettings const& fda_settings)
 : Base(result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES, syslen),
   result_type(result_type),
   syslen(syslen),
   distributed_forces(syslen, fda_settings),
   fda_settings(fda_settings)
{
    result_file << std::scientific << std::setprecision(6);
    if (PF_or_PS_mode()) make_backup(result_filename.c_str());
    if (fda_settings.binary_result_file) {
    	result_file.open(result_filename, std::ifstream::binary);
    	char b = 'b';
		result_file.write(&b, 1);
	    if (result_type == ResultType::PUNCTUAL_STRESS) {
	    	result_file.write(reinterpret_cast<char*>(&syslen), sizeof(uint));
	    }
    } else {
    	result_file.open(result_filename);
    	result_file << result_type << std::endl;
    }
    write_compat_header(1);
}

template <class Base>
void FDABase<Base>::write_frame(gmx::PaddedHostVector<gmx::RVec> const& x, const matrix box, int nsteps)
{
    switch (fda_settings.one_pair) {
        case OnePair::DETAILED:
            switch (result_type) {
                case ResultType::NO:
                    // do nothing
                    break;
                case ResultType::PAIRWISE_FORCES_VECTOR:
                    write_frame_detailed(x, box, true, nsteps);
                    break;
                case ResultType::PAIRWISE_FORCES_SCALAR:
                    write_frame_detailed(x, box, false, nsteps);
                    break;
                case ResultType::PUNCTUAL_STRESS:
                    gmx_fatal(FARGS, "Punctual stress is not supported for detailed output.\n");
                    break;
                case ResultType::VIRIAL_STRESS:
                    gmx_fatal(FARGS, "Virial stress is not supported for detailed output.\n");
                    break;
                case ResultType::VIRIAL_STRESS_VON_MISES:
                    gmx_fatal(FARGS, "Virial stress von Mises is not supported for detailed output.\n");
                    break;
                case ResultType::COMPAT_BIN:
                    gmx_fatal(FARGS, "Compatibility binary mode is not supported for detailed output.\n");
                    break;
                case ResultType::COMPAT_ASCII:
                    gmx_fatal(FARGS, "Compatibility ascii mode is not supported for detailed output.\n");
                    break;
            }
            break;
        case OnePair::SUMMED:
            switch (result_type) {
                case ResultType::NO:
                    // do nothing
                    break;
                case ResultType::PAIRWISE_FORCES_VECTOR:
                    write_frame_summed(x, box, true, nsteps);
                    break;
                case ResultType::PAIRWISE_FORCES_SCALAR:
                    write_frame_summed(x, box, false, nsteps);
                    break;
                case ResultType::PUNCTUAL_STRESS:
                    write_total_forces(x);
                    break;
                case ResultType::VIRIAL_STRESS:
                    write_virial_sum();
                    break;
                case ResultType::VIRIAL_STRESS_VON_MISES:
                    write_virial_sum_von_mises();
                    break;
                case ResultType::COMPAT_BIN:
                    gmx_fatal(FARGS, "Compatibility binary mode is not supported for summed output.\n");
                    break;
                case ResultType::COMPAT_ASCII:
                    write_frame_summed_compat(x, box, nsteps);
                    break;
            }
            break;
    }
}

template <class Base>
void FDABase<Base>::write_frame_detailed(gmx::PaddedHostVector<gmx::RVec> const& x, const matrix box, bool print_vector, int nsteps)
{
	write_frame_number(nsteps);
    if (print_vector)
        distributed_forces.write_detailed_vector(result_file);
    else
        distributed_forces.write_detailed_scalar(result_file, x, box);
}

template <class Base>
void FDABase<Base>::write_frame_summed(gmx::PaddedHostVector<gmx::RVec> const& x, const matrix box, bool print_vector, int nsteps)
{
	write_frame_number(nsteps);
    if (print_vector)
        distributed_forces.write_summed_vector(result_file);
    else
        distributed_forces.write_summed_scalar(result_file, x, box);
}

template <class Base>
void FDABase<Base>::write_frame_scalar(int nsteps)
{
	write_frame_number(nsteps);
    distributed_forces.write_scalar(result_file);
}

template <>
void FDABase<Atom>::write_total_forces(gmx::PaddedHostVector<gmx::RVec> const& x)
{
    distributed_forces.write_total_forces(result_file, x);
}

template <>
void FDABase<Residue>::write_total_forces(gmx::PaddedHostVector<gmx::RVec> const& x)
{
    distributed_forces.write_total_forces(result_file, x, fda_settings.normalize_psr);
}

template <class Base>
void FDABase<Base>::write_compat_header(int nsteps)
{
    if (!PF_or_PS_mode() or !compatibility_mode()) return;

    result_file << "<begin_block>" << std::endl;
    result_file << "; Forcemat version " << FDASettings::compat_fm_version << std::endl;
    result_file << "; Matrix containing pairwise forces." << std::endl;
    result_file << "; Matrix dimension " << syslen << " x " << syslen << std::endl;
    result_file << "version=" << FDASettings::compat_fm_version << std::endl;
    result_file << "groupname=" << fda_settings.groupname << std::endl;
    result_file << "writefreq=1" << std::endl;
    // Reserve place for up to 8 digits for nsteps
    result_file << "nsteps=" << std::setfill('0') << std::setw(8) << nsteps << std::endl;
    result_file << "sysanr=" << syslen << std::endl;
    result_file << "fmdim=" << syslen << std::endl;
    result_file << "intsize=" << sizeof(int) << std::endl;
    result_file << "realsize=" << sizeof(real) << std::endl;
    result_file << "<end_block>" << std::endl;
}

template <class Base>
void FDABase<Base>::write_frame_scalar_compat(int nsteps)
{
    if (result_type == ResultType::COMPAT_ASCII) {
        result_file << "<begin_block>" << std::endl;
        result_file << nsteps << std::endl;
        distributed_forces.write_scalar_compat_ascii(result_file);
        result_file << "<end_block>" << std::endl;
    } else {
        result_file.write(reinterpret_cast<const char *>(&nsteps), sizeof(nsteps));
        distributed_forces.write_scalar_compat_bin(result_file);
    }
}

template <class Base>
void FDABase<Base>::write_frame_summed_compat(gmx::PaddedHostVector<gmx::RVec> const& x, const matrix box, int nsteps)
{
    if (result_type == ResultType::COMPAT_ASCII) {
        result_file << "<begin_block>" << std::endl;
        result_file << nsteps << std::endl;
        distributed_forces.write_summed_compat_ascii(result_file, x, box);
        result_file << "<end_block>" << std::endl;
    } else {
        result_file.write(reinterpret_cast<const char *>(&nsteps), sizeof(nsteps));
        distributed_forces.write_summed_compat_bin(result_file, x, box);
    }
}

template <>
void FDABase<Atom>::write_virial_sum()
{
    bool first = true;
    for (auto const& v : virial_stress) {
        if (!first) result_file << " ";
        else first = false;
        result_file << -v(XX, XX) << " " << -v(YY, YY) << " " << -v(ZZ, ZZ) << " "
                    << -v(XX, YY) << " " << -v(XX, ZZ) << " " << -v(YY, ZZ);
    }
    result_file << std::endl;
}

template <>
void FDABase<Residue>::write_virial_sum()
{}

template <>
void FDABase<Atom>::write_virial_sum_von_mises()
{
    bool first = true;
    for (auto const& v : virial_stress) {
        if (!first) result_file << " ";
        else first = false;
        result_file << tensor_to_vonmises(v);
    }
    result_file << std::endl;
}

template <>
void FDABase<Residue>::write_virial_sum_von_mises()
{}

template <class Base>
void FDABase<Base>::write_frame_number(int nsteps)
{
	if (!fda_settings.binary_result_file) result_file << "frame " << nsteps << std::endl;
}

/// template instantiation
template class FDABase<Atom>;
template class FDABase<Residue>;

} // namespace fda
