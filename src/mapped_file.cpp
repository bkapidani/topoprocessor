/*
 * This source file is part of EMT, the ElectroMagneticTool.
 *
 * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <cstring>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "mapped_file.h"

mapped_file::mapped_file()
    : _start(nullptr), _end(nullptr), _is_open(false)
{}

mapped_file::mapped_file(const std::string& name)
    : mapped_file()
{
    open(name);
}

mapped_file::~mapped_file()
{
    close();
}

void
mapped_file::open(const std::string& name)
{
    close();
    _is_open = load(name);
}

void
mapped_file::close() noexcept
{
    reset();
}

bool
mapped_file::load(const std::string& name)
{
    std::ifstream input(name, std::ios::binary | std::ios::ate);
    if (!input) {
        return false;
    }

    const std::streamoff length = input.tellg();
    if (length < 0 ||
        length > std::numeric_limits<std::streamsize>::max()) {
        return false;
    }

    const std::uintmax_t unsigned_length =
        static_cast<std::uintmax_t>(length);
    if (unsigned_length >= _buffer.max_size()) {
        return false;
    }
    const std::size_t size = static_cast<std::size_t>(unsigned_length);

    _buffer.assign(size + 1U, '\0');
    input.seekg(0, std::ios::beg);
    if (size != 0U &&
        !input.read(_buffer.data(), static_cast<std::streamsize>(size))) {
        reset();
        return false;
    }

    _start = _buffer.data();
    _end = _start + size;
    return true;
}

void
mapped_file::reset() noexcept
{
    std::vector<char>().swap(_buffer);
    _start = nullptr;
    _end = nullptr;
    _is_open = false;
}

bool
mapped_file::is_open() const noexcept
{
    return _is_open;
}

bool
mapped_file::end() const noexcept
{
    return !_is_open || _start == _end;
}

std::string
mapped_file::get_line()
{
    if (end()) {
        throw std::out_of_range("no more lines in mapped_file");
    }

    const char* const line_start = _start;
    const void* const newline =
        std::memchr(_start, '\n', static_cast<std::size_t>(_end - _start));
    if (newline == nullptr) {
        _start = _end;
        return std::string(line_start, _end);
    }

    const char* const line_end = static_cast<const char*>(newline);
    _start = line_end + 1;
    return std::string(line_start, line_end);
}

const char*
mapped_file::mem() const noexcept
{
    return _start;
}

