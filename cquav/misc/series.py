# Copyright 2022 James Adams
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import cadquery as cq

def series(shape, size = 5, length_offset=None, width_offset=None, height_offset=None, skip_last=0, skip_first=0, operation=None):
    series = cq.Assembly()
    length, width, height = __resolve_hit_box(shape)

    bounding_box = {
        "length":length,
        "width":width,
        "height":height
    }

    for i in range(size):
        #print(f'loop iteration {length} {width} {height}')
        x_coord = 0
        y_coord = 0
        z_coord = 0

        if length_offset != None:
            x_coord = i * (length + length_offset)
            if i != 0:
                bounding_box['length'] += length + length_offset

        if width_offset != None:
            y_coord = i * (width + width_offset)
            if i != 0:
                bounding_box['width'] += width + width_offset

        if height_offset != None:
            z_coord = i * (height + height_offset)
            if i != 0:
                bounding_box['height'] += height + height_offset

        if i < size-skip_last:
            if i >= skip_first:
                series_shape = shape
                if operation:
                    series_shape = operation(shape, size, i, bounding_box)
                series.add(series_shape,  loc=cq.Location(cq.Vector(x_coord, y_coord, z_coord)))
            else:
                print('skipping first series tile')
        else:
            print('skipping last series tile')

    comp = series.toCompound()
    work = cq.Workplane("XZ").center(0, 0).workplane()
    work.add(comp)

    #print(f'I think the bounding box is', bounding_box)
    # zero out the offset caused by the first node
    work = work.translate((length/2,width/2,height/2))

    # center based on bounding box
    work = work.translate((-1*(bounding_box['length']/2), -1*(bounding_box['width']/2), -1*(bounding_box['height']/2)))
    meta = {'type':'series', 'length':bounding_box['length'], 'width':bounding_box['width'], 'height':bounding_box['height']}
    work.metadata = meta
    return work

def __resolve_hit_box(shape):
    attributes = dir(shape)
    if 'metadata' in attributes:
        meta = shape.metadata
        length = meta['length']
        width = meta['width']
        height = meta['height']
    else:
        bounds = shape.val().BoundingBox()
        if bounds:
            length = bounds.xlen
            width = bounds.ylen
            height = bounds.zlen
        else:
            raise Exception('Could not resolve shape metadata for series. OR bounding box')
    return length, width, height
