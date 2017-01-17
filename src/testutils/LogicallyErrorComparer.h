/*
 * LogicallyErrorComparer.h
 *
 *  Created on: Oct 9, 2014
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef LOGICALLYERRORCOMPARER_H_
#define LOGICALLYERRORCOMPARER_H_

struct LogicallyEqualComparerBase
{
	LogicallyEqualComparerBase(double error_factor) : error_factor(error_factor) {}

	double error_factor;
};

template <bool weight_by_magnitude = true, bool ignore_sign = false>
struct LogicallyEqualComparer : public LogicallyEqualComparerBase
{
	LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

	template <class T>
	bool operator () (T a, T b) const
	{
		return a==b || std::abs(a-b) < std::abs(std::min(a,b)) * std::numeric_limits<real>::epsilon() * error_factor;
	}
};

template <>
struct LogicallyEqualComparer<true,true> : public LogicallyEqualComparerBase
{
	LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

	template <class T>
	bool operator () (T a, T b) const
	{
		return a==b || std::abs(a) - std::abs(b) < std::abs(std::min(a,b)) * std::numeric_limits<real>::epsilon() * error_factor;
	}
};

template <>
struct LogicallyEqualComparer<false,false> : public LogicallyEqualComparerBase
{
	LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

	template <class T>
	bool operator () (T a, T b) const
	{
		return a==b || std::abs(a-b) < std::numeric_limits<real>::epsilon() * error_factor;
	}
};

template <>
struct LogicallyEqualComparer<false,true> : public LogicallyEqualComparerBase
{
	LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

	template <class T>
	bool operator () (T a, T b) const
	{
		return a==b || std::abs(a) - std::abs(b) < std::numeric_limits<real>::epsilon() * error_factor;
	}
};

#endif /* LOGICALLYERRORCOMPARER_H_ */
