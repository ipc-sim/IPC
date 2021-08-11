//
//  Triplet.hpp
//  IPC
//
//  Created by Minchen Li on 9/30/18.
//

#ifndef Triplet_h
#define Triplet_h

#include <cassert>

#include <iostream>

#include <boost/functional/hash.hpp>

namespace IPC {

class Triplet {
public:
    std::array<int, 3> key;

    Triplet(const int* p_key)
    {
        key[0] = p_key[0];
        key[1] = p_key[1];
        key[2] = p_key[2];
    }
    Triplet(int key0, int key1, int key2)
    {
        key[0] = key0;
        key[1] = key1;
        key[2] = key2;
    }

    int operator[](int i) const
    {
        assert(0 <= i && i <= 2);
        return key[i];
    }

    bool operator<(const Triplet& right) const
    {
        if (key[0] < right.key[0]) {
            return true;
        }
        else if (key[0] == right.key[0]) {
            if (key[1] < right.key[1]) {
                return true;
            }
            else if (key[1] == right.key[1]) {
                if (key[2] < right.key[2]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator==(const Triplet& right) const
    {
        return key[0] == right[0] && key[1] == right[1] && key[2] == right[2];
    }
};

} // namespace IPC

namespace std {
template <>
struct hash<IPC::Triplet> {
    size_t operator()(IPC::Triplet const& triplet) const noexcept
    {
        return boost::hash_range(triplet.key.begin(), triplet.key.end());
    }
};
}

inline std::ostream& operator<<(std::ostream& os, const IPC::Triplet& t)
{
    os << t.key[0] << ' ' << t.key[1] << ' ' << t.key[2];
    return os;
}

#endif /* Triplet_h */
