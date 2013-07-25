#include "common/ConfigMap.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <gtest/gtest.h>

#include <sstream>
#include <string>

namespace ba = boost::archive;
namespace bs = boost::serialization;

namespace {
    template<typename ArchOut, typename ArchIn>
    struct Librarian {
        typedef ArchOut output_type;
        typedef ArchIn input_type;
    };

    typedef ConfigMap<std::string, int>::type MapType;
}

template<typename LibrarianType>
class TestConfigMap : public ::testing::Test {
public:

    void SetUp() {
        _m["one"] = 1;
        _m["two"] = 2;
        _m["three"] = 3;
    }

    // Archivers get mad if they get the map as a non-const ref
    MapType const& testMap() const {
        return _m;
    }

private:
    MapType _m;
};

typedef ::testing::Types<
    Librarian<ba::xml_oarchive, ba::xml_iarchive>,
    Librarian<ba::text_oarchive, ba::text_iarchive>,
    Librarian<ba::binary_oarchive, ba::binary_iarchive>
    > MyTypes;

TYPED_TEST_CASE(TestConfigMap, MyTypes);

TYPED_TEST(TestConfigMap, serialize) {
    std::stringstream ss;

    typename TypeParam::output_type out(ss);
    out << bs::make_nvp("map", this->testMap());

    MapType mCopy;
    typename TypeParam::input_type in(ss);
    in >> bs::make_nvp("map", mCopy);

    EXPECT_EQ(this->testMap(), mCopy);
}
