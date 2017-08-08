/*
 * =====================================================================================
 *
 *       Filename:  namespace_jj.cpp
 *
 *    Description:  implementation of functions in new namespace
 *
 *        Version:  1.0
 *        Created:  07/28/2015 09:30:23 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Chai (jjchai), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#include "namespace_jj.h"
namespace jj{
	template <typename T>
	T complement(const T & i)
	{
		return (-i);
	}
	Edge * complement( Edge * edge)
	{
		return edge->getReverseEdge();
	}
}
