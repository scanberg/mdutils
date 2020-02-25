#pragma once

struct TrajectoryFrame;

struct TrajectoryLoader {
    virtual bool load_frame(TrajectoryFrame* frame) = 0;
};