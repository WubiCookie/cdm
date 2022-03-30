/* cdm_concepts - v0.1.0 - C++20 utility concepts - https://github.com/WubiCookie/cdm
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

#ifndef CDM_CONCEPTS_HPP
#define CDM_CONCEPTS_HPP

#include <concepts>
#include <type_traits>

namespace cdm
{
template <typename T>
concept resizable = requires(T& t)
{
	t.resize(0);
};

template <typename T>
concept trivial = std::is_trivial_v<T>;

template <typename T>
concept standard_layout = std::is_standard_layout_v<T>;

template <typename T>
concept arithmetic = std::integral<T> || std::floating_point<T>;
}  // namespace cdm

#endif  // CDM_CONCEPTS_HPP
