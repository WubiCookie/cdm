/* cdm_utils v2.0.0
   C++20 utility library
   https://github.com/WubiCookie/cdm
   no warranty implied; use at your own risk

LICENSE

       DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004

Copyright (C) 2022 Charles Seizilles de Mazancourt <charles DOT de DOT
mazancourt AT hotmail DOT fr>

Everyone is permitted to copy and distribute verbatim or modified
copies of this license document, and changing it is allowed as long
as the name is changed.

           DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
  TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

 0. You just DO WHAT THE FUCK YOU WANT TO.

CREDITS

Written by Charles Seizilles de Mazancourt
*/

#ifndef CDM_UTILS_HPP
#define CDM_UTILS_HPP 1

#include <cdm_concepts.hpp>

namespace cdm
{
template <arithmetic T>
class counter
{
	T m_begin;
	T m_end;
	T m_current = m_begin;

public:
	counter(T begin, T end) : m_begin{begin}, m_end{end} {}
	counter(T count) : counter{T(0), count} {}

	class iterator
	{
		friend counter;

		T m_current = T(0);
		T m_end;

		iterator(T begin, T end) : m_current{begin}, m_end{end} {}

	public:
		iterator(const iterator&) = default;
		iterator(iterator&&) = default;
		iterator& operator=(const iterator&) = default;
		iterator& operator=(iterator&&) = default;

		iterator& operator++()
		{
			++m_current;
			return *this;
		}
		iterator operator++(int)
		{
			auto copy = *this;
			++m_current;
			return copy;
		}

		bool operator==(iterator o) const { return m_current == o.m_current; }
		bool operator!=(iterator o) const { return not operator==(o); }

		T& operator*() { return m_current; }
		const T& operator*() const { return m_current; }
	};

	iterator begin() const { return iterator{m_begin, m_end}; }
	iterator end() const { return iterator{m_end, m_end}; }
	iterator cbegin() const { return iterator{m_begin, m_end}; }
	iterator cend() const { return iterator{m_end, m_end}; }
};
}  // namespace cdm

#endif  // CDM_UTILS_HPP
