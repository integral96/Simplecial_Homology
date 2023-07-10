#ifndef GROUP_COMM_HPP
#define GROUP_COMM_HPP

#include <base_math.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/difference.hpp>
#include <boost/hana/mult.hpp>

namespace hana = boost::hana;
namespace mpl  = boost::mpl;
template<char C> requires (C == 'e' || C == 'i' || C == 'j' || C == 'k' || C == 'l')
struct quaternion {};
template<>
struct quaternion<'e'> {
    static constexpr char ord = 'e';
    static constexpr auto value = hana::make_tuple(hana::int_c<1>, hana::int_c<1>, hana::int_c<1>, hana::int_c<1>);
};
template<>
struct quaternion<'i'> {
    static constexpr char ord = 'i';
    static constexpr auto value = hana::make_tuple(hana::int_c<1>, hana::int_c<0>, hana::int_c<0>, hana::int_c<0>);
};
template<>
struct quaternion<'j'> {
    static constexpr char ord = 'j';
    static constexpr auto value = hana::make_tuple(hana::int_c<0>, hana::int_c<1>, hana::int_c<0>, hana::int_c<0>);
};
template<>
struct quaternion<'k'> {
    static constexpr char ord = 'k';
    static constexpr auto value = hana::make_tuple(hana::int_c<0>, hana::int_c<0>, hana::int_c<1>, hana::int_c<0>);
};
template<>
struct quaternion<'l'> {
    static constexpr char ord = 'l';
    static constexpr auto value = hana::make_tuple(hana::int_c<0>, hana::int_c<0>, hana::int_c<0>, hana::int_c<1>);
};

template<class left, class right>
struct mult : mpl::if_c<left::ord == right::ord, quaternion<'e'>,
                        typename mpl::if_c<left::ord == 'i' && right::ord == 'j', quaternion<'k'>,
                        typename mpl::if_c<left::ord == 'e', quaternion<right::ord>, quaternion<left::ord>>::type
                        >::type

                        >::type {
};


template<typename ...Args>
inline constexpr auto make_group(Args&&... args) {
    return hana::make<hana::tuple_tag>(std::forward<Args>(args)...);
}
template<int N, typename Grp_type> requires (N == decltype(hana::size(std::declval<Grp_type>()))::value)
class group_comm {
    using value_type = Grp_type;
    Grp_type group_;
public:
    group_comm() {}
    group_comm(const Grp_type& grp) : group_(grp) {

    }

    friend std::ostream& operator << (std::ostream& out, const group_comm& cmpl) {
        out << "<";
        hana::for_each(cmpl.group_, [](const auto& x) {
//            static std::string name = boost::core::demangle(typeid(x).name());
            std::cout << x << " ";
        });
        out << ">";
        return out;
    }
};

#endif // GROUP_COMM_HPP
