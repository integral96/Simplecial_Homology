#ifndef GROUP_COMM_HPP
#define GROUP_COMM_HPP

#include <base_math.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/difference.hpp>

namespace hana = boost::hana;
namespace mpl  = boost::mpl;

template<int N, typename Field>
class group_comm {
    using value_type = Field;
    using group_type = boost::mp11::mp_repeat_c<hana::set<hana::type<Field>>, N>;
    group_type group_;
public:
    group_comm() {}
    template<typename ...Args, std::enable_if_t<(sizeof... (Args) == N) && is_same_v_<Field, remove_cvref_t<Args>...>>* = nullptr>
    group_comm(Args&&... args) : group_(hana::make_set(hana::make_type(std::forward<Args>(args))...)) {

    }
    template<int M, typename Field1> requires (std::is_same_v<Field1, Field> && M <= N)
    group_comm(const group_comm<M, Field1>& other) : group_(other.group_) {

    }
    template<int M, typename Field1> requires (std::is_same_v<Field1, Field> && M <= N)
    auto difference(const group_comm<M, Field1>& other) const {
        return hana::difference(group_, other.get_group());
    }
    const group_type& get_group() const {
        return group_;
    }
    friend std::ostream& operator << (std::ostream& out, const group_comm& cmpl) {
        out << "<";
        hana::for_each(cmpl.group_, [&](auto x) {
            static std::string name = boost::core::demangle(typeid(x).name());
            std::cout << name << std::endl;
        });
        out << ">";
        return out;
    }
};

#endif // GROUP_COMM_HPP
