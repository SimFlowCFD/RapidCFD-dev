#pragma once

#include "DeviceConfig.H"

namespace Foam {

class DeviceStream {
    cudaStream_t  stream_;

public:
    DeviceStream();
    ~DeviceStream();

    void synchronize() const;
    cudaStream_t operator()() const;
};

}

#include "DeviceStreamI.H"

