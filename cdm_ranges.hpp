/* cdm_ranges - v0.1.0 - ranges utility - https://github.com/WubiCookie/cdm
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

#ifndef CDM_RANGES_HPP
#define CDM_RANGES_HPP

#include <ranges>

namespace cdm
{
namespace ranges
{
template <std::ranges::contiguous_range R>
size_t raw_size(const R& r)
{
	return sizeof(*std::ranges::cdata(r)) * std::ranges::size(r);
}
}  // namespace ranges
}  // namespace cdm

#endif  // CDM_RANGES_HPP
