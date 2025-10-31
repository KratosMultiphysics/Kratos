/*
 * JSON schema validator for JSON for modern C++
 *
 * Copyright (c) 2016-2019 Patrick Boettcher <p@yai.se>.
 *
 * SPDX-License-Identifier: MIT
 *
 */
#ifndef NLOHMANN_JSON_SCHEMA_HPP__
#define NLOHMANN_JSON_SCHEMA_HPP__

#ifdef _WIN32
#	if defined(JSON_SCHEMA_VALIDATOR_EXPORTS)
#		define JSON_SCHEMA_VALIDATOR_API __declspec(dllexport)
#	elif defined(JSON_SCHEMA_VALIDATOR_IMPORTS)
#		define JSON_SCHEMA_VALIDATOR_API __declspec(dllimport)
#	else
#		define JSON_SCHEMA_VALIDATOR_API
#	endif
#else
#	define JSON_SCHEMA_VALIDATOR_API
#endif

#include <json/json.hpp>

#ifdef NLOHMANN_JSON_VERSION_MAJOR
#	if (NLOHMANN_JSON_VERSION_MAJOR * 10000 + NLOHMANN_JSON_VERSION_MINOR * 100 + NLOHMANN_JSON_VERSION_PATCH) < 30800
#		error "Please use this library with NLohmann's JSON version 3.8.0 or higher"
#	endif
#else
#	error "expected existing NLOHMANN_JSON_VERSION_MAJOR preproc variable, please update to NLohmann's JSON 3.8.0"
#endif

// make yourself a home - welcome to nlohmann's namespace
namespace nlohmann
{

// A class representing a JSON-URI for schemas derived from
// section 8 of JSON Schema: A Media Type for Describing JSON Documents
// draft-wright-json-schema-00
//
// New URIs can be derived from it using the derive()-method.
// This is useful for resolving refs or subschema-IDs in json-schemas.
//
// This is done implement the requirements described in section 8.2.
//
class JSON_SCHEMA_VALIDATOR_API json_uri
{
	std::string urn_;

	std::string scheme_;
	std::string authority_;
	std::string path_;

	json::json_pointer pointer_; // fragment part if JSON-Pointer
	std::string identifier_;     // fragment part if Locatation Independent ID

protected:
	// decodes a JSON uri and replaces all or part of the currently stored values
	void update(const std::string &uri);

	std::tuple<std::string, std::string, std::string, std::string, std::string> as_tuple() const
	{
		return std::make_tuple(urn_, scheme_, authority_, path_, identifier_ != "" ? identifier_ : pointer_.to_string());
	}

public:
	json_uri(const std::string &uri)
	{
		update(uri);
	}

	const std::string &scheme() const { return scheme_; }
	const std::string &authority() const { return authority_; }
	const std::string &path() const { return path_; }

	const json::json_pointer &pointer() const { return pointer_; }
	const std::string &identifier() const { return identifier_; }

	std::string fragment() const
	{
		if (identifier_ == "")
			return pointer_.to_string();
		else
			return identifier_;
	}

	std::string url() const { return location(); }
	std::string location() const;

	static std::string escape(const std::string &);

	// create a new json_uri based in this one and the given uri
	// resolves relative changes (pathes or pointers) and resets part if proto or hostname changes
	json_uri derive(const std::string &uri) const
	{
		json_uri u = *this;
		u.update(uri);
		return u;
	}

	// append a pointer-field to the pointer-part of this uri
	json_uri append(const std::string &field) const
	{
		if (identifier_ != "")
			return *this;

		json_uri u = *this;
		u.pointer_ /= field;
		return u;
	}

	std::string to_string() const;

	friend bool operator<(const json_uri &l, const json_uri &r)
	{
		return l.as_tuple() < r.as_tuple();
	}

	friend bool operator==(const json_uri &l, const json_uri &r)
	{
		return l.as_tuple() == r.as_tuple();
	}

	friend std::ostream &operator<<(std::ostream &os, const json_uri &u);
};

namespace json_schema
{

extern json draft7_schema_builtin;

typedef std::function<void(const json_uri & /*id*/, json & /*value*/)> schema_loader;
typedef std::function<void(const std::string & /*format*/, const std::string & /*value*/)> format_checker;
typedef std::function<void(const std::string & /*contentEncoding*/, const std::string & /*contentMediaType*/, const json & /*instance*/)> content_checker;

// Interface for validation error handlers
class JSON_SCHEMA_VALIDATOR_API error_handler
{
public:
	virtual ~error_handler() {}
	virtual void error(const json::json_pointer & /*ptr*/, const json & /*instance*/, const std::string & /*message*/) = 0;
};

class JSON_SCHEMA_VALIDATOR_API basic_error_handler : public error_handler
{
	bool error_{false};

public:
	void error(const json::json_pointer & /*ptr*/, const json & /*instance*/, const std::string & /*message*/) override
	{
		error_ = true;
	}

	virtual void reset() { error_ = false; }
	operator bool() const { return error_; }
};

/**
 * Checks validity of JSON schema built-in string format specifiers like 'date-time', 'ipv4', ...
 */
void JSON_SCHEMA_VALIDATOR_API default_string_format_check(const std::string &format, const std::string &value);

class root_schema;

class JSON_SCHEMA_VALIDATOR_API json_validator
{
	std::unique_ptr<root_schema> root_;

public:
	json_validator(schema_loader = nullptr, format_checker = nullptr, content_checker = nullptr);

	json_validator(const json &, schema_loader = nullptr, format_checker = nullptr, content_checker = nullptr);
	json_validator(json &&, schema_loader = nullptr, format_checker = nullptr, content_checker = nullptr);

	json_validator(json_validator &&);
	json_validator &operator=(json_validator &&);

	json_validator(json_validator const &) = delete;
	json_validator &operator=(json_validator const &) = delete;

	~json_validator();

	// insert and set the root-schema
	void set_root_schema(const json &);
	void set_root_schema(json &&);

	// validate a json-document based on the root-schema
	json validate(const json &) const;

	// validate a json-document based on the root-schema with a custom error-handler
	json validate(const json &, error_handler &, const json_uri &initial_uri = json_uri("#")) const;
};

} // namespace json_schema
} // namespace nlohmann

#endif /* NLOHMANN_JSON_SCHEMA_HPP__ */
#pragma once

#include <json/json.hpp>
#include <string>

namespace nlohmann
{
class JsonPatchFormatException : public std::exception
{
public:
	explicit JsonPatchFormatException(std::string msg)
	    : ex_{std::move(msg)} {}

	inline const char *what() const noexcept override final { return ex_.c_str(); }

private:
	std::string ex_;
};

class json_patch
{
public:
	json_patch() = default;
	json_patch(json &&patch);
	json_patch(const json &patch);

	json_patch &add(const json::json_pointer &, json value);
	json_patch &replace(const json::json_pointer &, json value);
	json_patch &remove(const json::json_pointer &);

	json &get_json() { return j_; }
	const json &get_json() const { return j_; }

	operator json() const { return j_; }

private:
	json j_ = nlohmann::json::array();

	static void validateJsonPatch(json const &patch);
};
} // namespace nlohmann
// #include "json-patch.hpp"

namespace
{

// originally from http://jsonpatch.com/, http://json.schemastore.org/json-patch
// with fixes
const nlohmann::json patch_schema = R"patch({
    "title": "JSON schema for JSONPatch files",
    "$schema": "http://json-schema.org/draft-04/schema#",
    "type": "array",

    "items": {
        "oneOf": [
            {
                "additionalProperties": false,
                "required": [ "value", "op", "path"],
                "properties": {
                    "path" : { "$ref": "#/definitions/path" },
                    "op": {
                        "description": "The operation to perform.",
                        "type": "string",
                        "enum": [ "add", "replace", "test" ]
                    },
                    "value": {
                        "description": "The value to add, replace or test."
                    }
                }
            },
            {
                "additionalProperties": false,
                "required": [ "op", "path"],
                "properties": {
                    "path" : { "$ref": "#/definitions/path" },
                    "op": {
                        "description": "The operation to perform.",
                        "type": "string",
                        "enum": [ "remove" ]
                    }
                }
            },
            {
                "additionalProperties": false,
                "required": [ "from", "op", "path" ],
                "properties": {
                    "path" : { "$ref": "#/definitions/path" },
                    "op": {
                        "description": "The operation to perform.",
                        "type": "string",
                        "enum": [ "move", "copy" ]
                    },
                    "from": {
                        "$ref": "#/definitions/path",
                        "description": "A JSON Pointer path pointing to the location to move/copy from."
                    }
                }
            }
        ]
    },
    "definitions": {
        "path": {
            "description": "A JSON Pointer path.",
            "type": "string"
        }
    }
})patch"_json;
} // namespace

namespace nlohmann
{

json_patch::json_patch(json &&patch)
    : j_(std::move(patch))
{
	validateJsonPatch(j_);
}

json_patch::json_patch(const json &patch)
    : j_(std::move(patch))
{
	validateJsonPatch(j_);
}

json_patch &json_patch::add(const json::json_pointer &ptr, json value)
{
	j_.push_back(json{{"op", "add"}, {"path", ptr.to_string()}, {"value", std::move(value)}});
	return *this;
}

json_patch &json_patch::replace(const json::json_pointer &ptr, json value)
{
	j_.push_back(json{{"op", "replace"}, {"path", ptr.to_string()}, {"value", std::move(value)}});
	return *this;
}

json_patch &json_patch::remove(const json::json_pointer &ptr)
{
	j_.push_back(json{{"op", "remove"}, {"path", ptr.to_string()}});
	return *this;
}

void json_patch::validateJsonPatch(json const &patch)
{
	// static put here to have it created at the first usage of validateJsonPatch
	static nlohmann::json_schema::json_validator patch_validator(patch_schema);

	patch_validator.validate(patch);

	for (auto const &op : patch)
		json::json_pointer(op["path"].get<std::string>());
}

} // namespace nlohmann
/*
 * JSON schema validator for JSON for modern C++
 *
 * Copyright (c) 2016-2019 Patrick Boettcher <p@yai.se>.
 *
 * SPDX-License-Identifier: MIT
 *
 */
#include <json/json.hpp>

namespace nlohmann
{
namespace json_schema
{

json draft7_schema_builtin = R"( {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "http://json-schema.org/draft-07/schema#",
    "title": "Core schema meta-schema",
    "definitions": {
        "schemaArray": {
            "type": "array",
            "minItems": 1,
            "items": { "$ref": "#" }
        },
        "nonNegativeInteger": {
            "type": "integer",
            "minimum": 0
        },
        "nonNegativeIntegerDefault0": {
            "allOf": [
                { "$ref": "#/definitions/nonNegativeInteger" },
                { "default": 0 }
            ]
        },
        "simpleTypes": {
            "enum": [
                "array",
                "boolean",
                "integer",
                "null",
                "number",
                "object",
                "string"
            ]
        },
        "stringArray": {
            "type": "array",
            "items": { "type": "string" },
            "uniqueItems": true,
            "default": []
        }
    },
    "type": ["object", "boolean"],
    "properties": {
        "$id": {
            "type": "string",
            "format": "uri-reference"
        },
        "$schema": {
            "type": "string",
            "format": "uri"
        },
        "$ref": {
            "type": "string",
            "format": "uri-reference"
        },
        "$comment": {
            "type": "string"
        },
        "title": {
            "type": "string"
        },
        "description": {
            "type": "string"
        },
        "default": true,
        "readOnly": {
            "type": "boolean",
            "default": false
        },
        "examples": {
            "type": "array",
            "items": true
        },
        "multipleOf": {
            "type": "number",
            "exclusiveMinimum": 0
        },
        "maximum": {
            "type": "number"
        },
        "exclusiveMaximum": {
            "type": "number"
        },
        "minimum": {
            "type": "number"
        },
        "exclusiveMinimum": {
            "type": "number"
        },
        "maxLength": { "$ref": "#/definitions/nonNegativeInteger" },
        "minLength": { "$ref": "#/definitions/nonNegativeIntegerDefault0" },
        "pattern": {
            "type": "string",
            "format": "regex"
        },
        "additionalItems": { "$ref": "#" },
        "items": {
            "anyOf": [
                { "$ref": "#" },
                { "$ref": "#/definitions/schemaArray" }
            ],
            "default": true
        },
        "maxItems": { "$ref": "#/definitions/nonNegativeInteger" },
        "minItems": { "$ref": "#/definitions/nonNegativeIntegerDefault0" },
        "uniqueItems": {
            "type": "boolean",
            "default": false
        },
        "contains": { "$ref": "#" },
        "maxProperties": { "$ref": "#/definitions/nonNegativeInteger" },
        "minProperties": { "$ref": "#/definitions/nonNegativeIntegerDefault0" },
        "required": { "$ref": "#/definitions/stringArray" },
        "additionalProperties": { "$ref": "#" },
        "definitions": {
            "type": "object",
            "additionalProperties": { "$ref": "#" },
            "default": {}
        },
        "properties": {
            "type": "object",
            "additionalProperties": { "$ref": "#" },
            "default": {}
        },
        "patternProperties": {
            "type": "object",
            "additionalProperties": { "$ref": "#" },
            "propertyNames": { "format": "regex" },
            "default": {}
        },
        "dependencies": {
            "type": "object",
            "additionalProperties": {
                "anyOf": [
                    { "$ref": "#" },
                    { "$ref": "#/definitions/stringArray" }
                ]
            }
        },
        "propertyNames": { "$ref": "#" },
        "const": true,
        "enum": {
            "type": "array",
            "items": true,
            "minItems": 1,
            "uniqueItems": true
        },
        "type": {
            "anyOf": [
                { "$ref": "#/definitions/simpleTypes" },
                {
                    "type": "array",
                    "items": { "$ref": "#/definitions/simpleTypes" },
                    "minItems": 1,
                    "uniqueItems": true
                }
            ]
        },
        "format": { "type": "string" },
        "contentMediaType": { "type": "string" },
        "contentEncoding": { "type": "string" },
        "if": { "$ref": "#" },
        "then": { "$ref": "#" },
        "else": { "$ref": "#" },
        "allOf": { "$ref": "#/definitions/schemaArray" },
        "anyOf": { "$ref": "#/definitions/schemaArray" },
        "oneOf": { "$ref": "#/definitions/schemaArray" },
        "not": { "$ref": "#" }
    },
    "default": true
} )"_json;
}
} // namespace nlohmann
/*
 * JSON schema validator for JSON for modern C++
 *
 * Copyright (c) 2016-2019 Patrick Boettcher <p@yai.se>.
 *
 * SPDX-License-Identifier: MIT
 *
 */
// #include <nlohmann/json-schema.hpp>


#include <sstream>

namespace nlohmann
{

void json_uri::update(const std::string &uri)
{
	std::string pointer = ""; // default pointer is document-root

	// first split the URI into location and pointer
	auto pointer_separator = uri.find('#');
	if (pointer_separator != std::string::npos) {    // and extract the pointer-string if found
		pointer = uri.substr(pointer_separator + 1); // remove #

		// unescape %-values IOW, decode JSON-URI-formatted JSON-pointer
		std::size_t pos = pointer.size() - 1;
		do {
			pos = pointer.rfind('%', pos);
			if (pos == std::string::npos)
				break;

			if (pos >= pointer.size() - 2) {
				pos--;
				continue;
			}

			std::string hex = pointer.substr(pos + 1, 2);
			char ascii = static_cast<char>(std::strtoul(hex.c_str(), nullptr, 16));
			pointer.replace(pos, 3, 1, ascii);

			pos--;
		} while (1);
	}

	auto location = uri.substr(0, pointer_separator);

	if (location.size()) { // a location part has been found

		// if it is an URN take it as it is
		if (location.find("urn:") == 0) {
			urn_ = location;

			// and clear URL members
			scheme_ = "";
			authority_ = "";
			path_ = "";

		} else { // it is an URL

			// split URL in protocol, hostname and path
			std::size_t pos = 0;
			auto proto = location.find("://", pos);
			if (proto != std::string::npos) { // extract the protocol

				urn_ = ""; // clear URN-member if URL is parsed

				scheme_ = location.substr(pos, proto - pos);
				pos = 3 + proto; // 3 == "://"

				auto authority = location.find("/", pos);
				if (authority != std::string::npos) { // and the hostname (no proto without hostname)
					authority_ = location.substr(pos, authority - pos);
					pos = authority;
				}
			}

			auto path = location.substr(pos);

			// URNs cannot of have paths
			if (urn_.size() && path.size())
				throw std::invalid_argument("Cannot add a path (" + path + ") to an URN URI (" + urn_ + ")");

			if (path[0] == '/') // if it starts with a / it is root-path
				path_ = path;
			else if (pos == 0) { // the URL contained only a path and the current path has no / at the end, strip last element until / and append
				auto last_slash = path_.rfind('/');
				path_ = path_.substr(0, last_slash) + '/' + path;
			} else // otherwise it is a subfolder
				path_.append(path);
		}
	}

	pointer_ = ""_json_pointer;
	identifier_ = "";

	if (pointer[0] == '/')
		pointer_ = json::json_pointer(pointer);
	else
		identifier_ = pointer;
}

std::string json_uri::location() const
{
	if (urn_.size())
		return urn_;

	std::stringstream s;

	if (scheme_.size() > 0)
		s << scheme_ << "://";

	s << authority_
	  << path_;

	return s.str();
}

std::string json_uri::to_string() const
{
	std::stringstream s;

	s << location() << " # ";

	if (identifier_ == "")
		s << pointer_.to_string();
	else
		s << identifier_;

	return s.str();
}

std::ostream &operator<<(std::ostream &os, const json_uri &u)
{
	return os << u.to_string();
}

std::string json_uri::escape(const std::string &src)
{
	std::vector<std::pair<std::string, std::string>> chars = {
	    {"~", "~0"},
	    {"/", "~1"}};

	std::string l = src;

	for (const auto &c : chars) {
		std::size_t pos = 0;
		do {
			pos = l.find(c.first, pos);
			if (pos == std::string::npos)
				break;
			l.replace(pos, 1, c.second);
			pos += c.second.size();
		} while (1);
	}

	return l;
}

} // namespace nlohmann
/*
 * JSON schema validator for JSON for modern C++
 *
 * Copyright (c) 2016-2019 Patrick Boettcher <p@yai.se>.
 *
 * SPDX-License-Identifier: MIT
 *
 */
// #include <nlohmann/json-schema.hpp>


// #include "json-patch.hpp"


#include <deque>
#include <memory>
#include <set>
#include <sstream>
#include <string>

using nlohmann::json;
using nlohmann::json_patch;
using nlohmann::json_uri;
using nlohmann::json_schema::root_schema;
using namespace nlohmann::json_schema;

#ifdef JSON_SCHEMA_BOOST_REGEX
#	include <boost/regex.hpp>
#	define REGEX_NAMESPACE boost
#elif defined(JSON_SCHEMA_NO_REGEX)
#	define NO_STD_REGEX
#else
#	include <regex>
#	define REGEX_NAMESPACE std
#endif

namespace
{

class schema
{
protected:
	root_schema *root_;
	json default_value_ = nullptr;

protected:
	virtual std::shared_ptr<schema> make_for_default_(
	    std::shared_ptr<::schema> & /* sch */,
	    root_schema * /* root */,
	    std::vector<nlohmann::json_uri> & /* uris */,
	    nlohmann::json & /* default_value */) const
	{
		return nullptr;
	};

public:
	virtual ~schema() = default;

	schema(root_schema *root)
	    : root_(root) {}

	virtual void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const = 0;

	virtual const json &default_value(const json::json_pointer &, const json &, error_handler &) const
	{
		return default_value_;
	}

	void set_default_value(const json &v) { default_value_ = v; }

	static std::shared_ptr<schema> make(json &schema,
	                                    root_schema *root,
	                                    const std::vector<std::string> &key,
	                                    std::vector<nlohmann::json_uri> uris);
};

class schema_ref : public schema
{
	const std::string id_;
	std::weak_ptr<schema> target_;
	std::shared_ptr<schema> target_strong_; // for references to references keep also the shared_ptr because
	                                        // no one else might use it after resolving

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const final
	{
		auto target = target_.lock();

		if (target)
			target->validate(ptr, instance, patch, e);
		else
			e.error(ptr, instance, "unresolved or freed schema-reference " + id_);
	}

	const json &default_value(const json::json_pointer &ptr, const json &instance, error_handler &e) const override final
	{
		if (!default_value_.is_null())
			return default_value_;

		auto target = target_.lock();
		if (target)
			return target->default_value(ptr, instance, e);

		e.error(ptr, instance, "unresolved or freed schema-reference " + id_);

		return default_value_;
	}

protected:
	virtual std::shared_ptr<schema> make_for_default_(
	    std::shared_ptr<::schema> &sch,
	    root_schema *root,
	    std::vector<nlohmann::json_uri> &uris,
	    nlohmann::json &default_value) const override
	{
		// create a new reference schema using the original reference (which will be resolved later)
		// to store this overloaded default value #209
		auto result = std::make_shared<schema_ref>(uris[0].to_string(), root);
		result->set_target(sch, true);
		result->set_default_value(default_value);
		return result;
	};

public:
	schema_ref(const std::string &id, root_schema *root)
	    : schema(root), id_(id) {}

	const std::string &id() const { return id_; }

	void set_target(const std::shared_ptr<schema> &target, bool strong = false)
	{
		target_ = target;
		if (strong)
			target_strong_ = target;
	}
};

} // namespace

namespace nlohmann
{
namespace json_schema
{

class root_schema
{
	schema_loader loader_;
	format_checker format_check_;
	content_checker content_check_;

	std::shared_ptr<schema> root_;

	struct schema_file {
		std::map<std::string, std::shared_ptr<schema>> schemas;
		std::map<std::string, std::shared_ptr<schema_ref>> unresolved; // contains all unresolved references from any other file seen during parsing
		json unknown_keywords;
	};

	// location as key
	std::map<std::string, schema_file> files_;

	schema_file &get_or_create_file(const std::string &loc)
	{
		auto file = files_.lower_bound(loc);
		if (file != files_.end() && !(files_.key_comp()(loc, file->first)))
			return file->second;
		else
			return files_.insert(file, {loc, {}})->second;
	}

public:
	root_schema(schema_loader &&loader,
	            format_checker &&format,
	            content_checker &&content)

	    : loader_(std::move(loader)),
	      format_check_(std::move(format)),
	      content_check_(std::move(content))
	{
	}

	format_checker &format_check() { return format_check_; }
	content_checker &content_check() { return content_check_; }

	void insert(const json_uri &uri, const std::shared_ptr<schema> &s)
	{
		auto &file = get_or_create_file(uri.location());
		auto sch = file.schemas.lower_bound(uri.fragment());
		if (sch != file.schemas.end() && !(file.schemas.key_comp()(uri.fragment(), sch->first))) {
			throw std::invalid_argument("schema with " + uri.to_string() + " already inserted");
			return;
		}

		file.schemas.insert({uri.fragment(), s});

		// was someone referencing this newly inserted schema?
		auto unresolved = file.unresolved.find(uri.fragment());
		if (unresolved != file.unresolved.end()) {
			unresolved->second->set_target(s);
			file.unresolved.erase(unresolved);
		}
	}

	void insert_unknown_keyword(const json_uri &uri, const std::string &key, json &value)
	{
		auto &file = get_or_create_file(uri.location());
		auto new_uri = uri.append(key);
		auto fragment = new_uri.pointer();

		// is there a reference looking for this unknown-keyword, which is thus no longer a unknown keyword but a schema
		auto unresolved = file.unresolved.find(fragment.to_string());
		if (unresolved != file.unresolved.end())
			schema::make(value, this, {}, {{new_uri}});
		else { // no, nothing ref'd it, keep for later

			// need to create an object for each reference-token in the
			// JSON-Pointer When not existing, a stringified integer reference
			// token (e.g. "123") in the middle of the pointer will be
			// interpreted a an array-index and an array will be created.

			// json_pointer's reference_tokens is private - get them
			std::deque<std::string> ref_tokens;
			auto uri_pointer = uri.pointer();
			while (!uri_pointer.empty()) {
				ref_tokens.push_front(uri_pointer.back());
				uri_pointer.pop_back();
			}

			// for each token create an object, if not already existing
			auto unk_kw = &file.unknown_keywords;
			for (auto &rt : ref_tokens) {
				// create a json_pointer from rt as rt can be an stringified integer doing find on an array won't work
				json::json_pointer rt_ptr{"/" + rt};
				if (unk_kw->contains(rt_ptr) == false)
					(*unk_kw)[rt] = json::object();
				unk_kw = &(*unk_kw)[rt_ptr];
			}
			(*unk_kw)[key] = value;
		}

		// recursively add possible subschemas of unknown keywords
		if (value.type() == json::value_t::object)
			for (auto &subsch : value.items())
				insert_unknown_keyword(new_uri, subsch.key(), subsch.value());
	}

	std::shared_ptr<schema> get_or_create_ref(const json_uri &uri)
	{
		auto &file = get_or_create_file(uri.location());

		// existing schema
		auto sch = file.schemas.find(uri.fragment());
		if (sch != file.schemas.end())
			return sch->second;

		// referencing an unknown keyword, turn it into schema
		//
		// an unknown keyword can only be referenced by a json-pointer,
		// not by a plain name fragment
		if (!uri.pointer().to_string().empty()) {
			bool contains_pointer = file.unknown_keywords.contains(uri.pointer());
			if (contains_pointer) {
				auto &subschema = file.unknown_keywords.at(uri.pointer());
				auto s = schema::make(subschema, this, {}, {{uri}});
				if (s) { // if schema is valid (non-null)
					file.unknown_keywords.erase(uri.fragment());
					return s;
				}
			}
		}

		// get or create a schema_ref
		auto r = file.unresolved.lower_bound(uri.fragment());
		if (r != file.unresolved.end() && !(file.unresolved.key_comp()(uri.fragment(), r->first))) {
			return r->second; // unresolved, already seen previously - use existing reference
		} else {
			return file.unresolved.insert(r,
			                              {uri.fragment(), std::make_shared<schema_ref>(uri.to_string(), this)})
			    ->second; // unresolved, create reference
		}
	}

	void set_root_schema(json sch)
	{
		files_.clear();
		root_ = schema::make(sch, this, {}, {{"#"}});

		// load all files which have not yet been loaded
		do {
			bool new_schema_loaded = false;

			// files_ is modified during parsing, iterators are invalidated
			std::vector<std::string> locations;
			for (auto &file : files_)
				locations.push_back(file.first);

			for (auto &loc : locations) {
				if (files_[loc].schemas.size() == 0) { // nothing has been loaded for this file
					if (loader_) {
						json loaded_schema;

						loader_(loc, loaded_schema);

						schema::make(loaded_schema, this, {}, {{loc}});
						new_schema_loaded = true;
					} else {
						throw std::invalid_argument("external schema reference '" + loc + "' needs loading, but no loader callback given");
					}
				}
			}

			if (!new_schema_loaded) // if no new schema loaded, no need to try again
				break;
		} while (1);

		for (const auto &file : files_) {
			if (file.second.unresolved.size() != 0) {
				// Build a representation of the undefined
				// references as a list of comma-separated strings.
				auto n_urefs = file.second.unresolved.size();
				std::string urefs = "[";

				decltype(n_urefs) counter = 0;
				for (const auto &p : file.second.unresolved) {
					urefs += p.first;

					if (counter != n_urefs - 1u) {
						urefs += ", ";
					}

					++counter;
				}

				urefs += "]";

				throw std::invalid_argument("after all files have been parsed, '" +
				                            (file.first == "" ? "<root>" : file.first) +
				                            "' has still the following undefined references: " + urefs);
			}
		}
	}

	void validate(const json::json_pointer &ptr,
	              const json &instance,
	              json_patch &patch,
	              error_handler &e,
	              const json_uri &initial) const
	{
		if (!root_) {
			e.error(ptr, "", "no root schema has yet been set for validating an instance");
			return;
		}

		auto file_entry = files_.find(initial.location());
		if (file_entry == files_.end()) {
			e.error(ptr, "", "no file found serving requested root-URI. " + initial.location());
			return;
		}

		auto &file = file_entry->second;
		auto sch = file.schemas.find(initial.fragment());
		if (sch == file.schemas.end()) {
			e.error(ptr, "", "no schema find for request initial URI: " + initial.to_string());
			return;
		}

		sch->second->validate(ptr, instance, patch, e);
	}
};

} // namespace json_schema
} // namespace nlohmann

namespace
{

class first_error_handler : public error_handler
{
public:
	bool error_{false};
	json::json_pointer ptr_;
	json instance_;
	std::string message_;

	void error(const json::json_pointer &ptr, const json &instance, const std::string &message) override
	{
		if (*this)
			return;
		error_ = true;
		ptr_ = ptr;
		instance_ = instance;
		message_ = message;
	}

	operator bool() const { return error_; }
};

class logical_not : public schema
{
	std::shared_ptr<schema> subschema_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const final
	{
		first_error_handler esub;
		subschema_->validate(ptr, instance, patch, esub);

		if (!esub)
			e.error(ptr, instance, "the subschema has succeeded, but it is required to not validate");
	}

	const json &default_value(const json::json_pointer &ptr, const json &instance, error_handler &e) const override
	{
		return subschema_->default_value(ptr, instance, e);
	}

public:
	logical_not(json &sch,
	            root_schema *root,
	            const std::vector<nlohmann::json_uri> &uris)
	    : schema(root)
	{
		subschema_ = schema::make(sch, root, {"not"}, uris);
	}
};

enum logical_combination_types {
	allOf,
	anyOf,
	oneOf
};

class logical_combination_error_handler : public error_handler
{
public:
	struct error_entry {
		json::json_pointer ptr_;
		json instance_;
		std::string message_;
	};

	std::vector<error_entry> error_entry_list_;

	void error(const json::json_pointer &ptr, const json &instance, const std::string &message) override
	{
		error_entry_list_.push_back(error_entry{ptr, instance, message});
	}

	void propagate(error_handler &e, const std::string &prefix) const
	{
		for (const error_entry &entry : error_entry_list_)
			e.error(entry.ptr_, entry.instance_, prefix + entry.message_);
	}

	operator bool() const { return !error_entry_list_.empty(); }
};

template <enum logical_combination_types combine_logic>
class logical_combination : public schema
{
	std::vector<std::shared_ptr<schema>> subschemata_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const final
	{
		size_t count = 0;
		logical_combination_error_handler error_summary;

		for (std::size_t index = 0; index < subschemata_.size(); ++index) {
			const std::shared_ptr<schema> &s = subschemata_[index];
			logical_combination_error_handler esub;
			auto oldPatchSize = patch.get_json().size();
			s->validate(ptr, instance, patch, esub);
			if (!esub)
				count++;
			else {
				patch.get_json().get_ref<nlohmann::json::array_t &>().resize(oldPatchSize);
				esub.propagate(error_summary, "case#" + std::to_string(index) + "] ");
			}

			if (is_validate_complete(instance, ptr, e, esub, count, index))
				return;
		}

		if (count == 0) {
			e.error(ptr, instance, "no subschema has succeeded, but one of them is required to validate. Type: " + key + ", number of failed subschemas: " + std::to_string(subschemata_.size()));
			error_summary.propagate(e, "[combination: " + key + " / ");
		}
	}

	// specialized for each of the logical_combination_types
	static const std::string key;
	static bool is_validate_complete(const json &, const json::json_pointer &, error_handler &, const logical_combination_error_handler &, size_t, size_t);

public:
	logical_combination(json &sch,
	                    root_schema *root,
	                    const std::vector<nlohmann::json_uri> &uris)
	    : schema(root)
	{
		size_t c = 0;
		for (auto &subschema : sch)
			subschemata_.push_back(schema::make(subschema, root, {key, std::to_string(c++)}, uris));

		// value of allOf, anyOf, and oneOf "MUST be a non-empty array"
		// TODO error/throw? when subschemata_.empty()
	}
};

template <>
const std::string logical_combination<allOf>::key = "allOf";
template <>
const std::string logical_combination<anyOf>::key = "anyOf";
template <>
const std::string logical_combination<oneOf>::key = "oneOf";

template <>
bool logical_combination<allOf>::is_validate_complete(const json &, const json::json_pointer &, error_handler &e, const logical_combination_error_handler &esub, size_t, size_t current_schema_index)
{
	if (esub) {
		e.error(esub.error_entry_list_.front().ptr_, esub.error_entry_list_.front().instance_, "at least one subschema has failed, but all of them are required to validate - " + esub.error_entry_list_.front().message_);
		esub.propagate(e, "[combination: allOf / case#" + std::to_string(current_schema_index) + "] ");
	}
	return esub;
}

template <>
bool logical_combination<anyOf>::is_validate_complete(const json &, const json::json_pointer &, error_handler &, const logical_combination_error_handler &, size_t count, size_t)
{
	return count == 1;
}

template <>
bool logical_combination<oneOf>::is_validate_complete(const json &instance, const json::json_pointer &ptr, error_handler &e, const logical_combination_error_handler &, size_t count, size_t)
{
	if (count > 1)
		e.error(ptr, instance, "more than one subschema has succeeded, but exactly one of them is required to validate");
	return count > 1;
}

class type_schema : public schema
{
	std::vector<std::shared_ptr<schema>> type_;
	std::pair<bool, json> enum_, const_;
	std::vector<std::shared_ptr<schema>> logic_;

	static std::shared_ptr<schema> make(json &schema,
	                                    json::value_t type,
	                                    root_schema *,
	                                    const std::vector<nlohmann::json_uri> &,
	                                    std::set<std::string> &);

	std::shared_ptr<schema> if_, then_, else_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const override final
	{
		// depending on the type of instance run the type specific validator - if present
		auto type = type_[static_cast<uint8_t>(instance.type())];

		if (type)
			type->validate(ptr, instance, patch, e);
		else
			e.error(ptr, instance, "unexpected instance type");

		if (enum_.first) {
			bool seen_in_enum = false;
			for (auto &v : enum_.second)
				if (instance == v) {
					seen_in_enum = true;
					break;
				}

			if (!seen_in_enum)
				e.error(ptr, instance, "instance not found in required enum");
		}

		if (const_.first &&
		    const_.second != instance)
			e.error(ptr, instance, "instance not const");

		for (auto l : logic_)
			l->validate(ptr, instance, patch, e);

		if (if_) {
			first_error_handler err;

			if_->validate(ptr, instance, patch, err);
			if (!err) {
				if (then_)
					then_->validate(ptr, instance, patch, e);
			} else {
				if (else_)
					else_->validate(ptr, instance, patch, e);
			}
		}
		if (instance.is_null()) {
			patch.add(nlohmann::json::json_pointer{}, default_value_);
		}
	}

protected:
	virtual std::shared_ptr<schema> make_for_default_(
	    std::shared_ptr<::schema> & /* sch */,
	    root_schema * /* root */,
	    std::vector<nlohmann::json_uri> & /* uris */,
	    nlohmann::json &default_value) const override
	{
		auto result = std::make_shared<type_schema>(*this);
		result->set_default_value(default_value);
		return result;
	};

public:
	type_schema(json &sch,
	            root_schema *root,
	            const std::vector<nlohmann::json_uri> &uris)
	    : schema(root), type_(static_cast<uint8_t>(json::value_t::discarded) + 1)
	{
		// association between JSON-schema-type and NLohmann-types
		static const std::vector<std::pair<std::string, json::value_t>> schema_types = {
		    {"null", json::value_t::null},
		    {"object", json::value_t::object},
		    {"array", json::value_t::array},
		    {"string", json::value_t::string},
		    {"boolean", json::value_t::boolean},
		    {"integer", json::value_t::number_integer},
		    {"number", json::value_t::number_float},
		};

		std::set<std::string> known_keywords;

		auto attr = sch.find("type");
		if (attr == sch.end()) // no type field means all sub-types possible
			for (auto &t : schema_types)
				type_[static_cast<uint8_t>(t.second)] = type_schema::make(sch, t.second, root, uris, known_keywords);
		else {
			switch (attr.value().type()) { // "type": "type"

			case json::value_t::string: {
				auto schema_type = attr.value().get<std::string>();
				for (auto &t : schema_types)
					if (t.first == schema_type)
						type_[static_cast<uint8_t>(t.second)] = type_schema::make(sch, t.second, root, uris, known_keywords);
			} break;

			case json::value_t::array: // "type": ["type1", "type2"]
				for (auto &array_value : attr.value()) {
					auto schema_type = array_value.get<std::string>();
					for (auto &t : schema_types)
						if (t.first == schema_type)
							type_[static_cast<uint8_t>(t.second)] = type_schema::make(sch, t.second, root, uris, known_keywords);
				}
				break;

			default:
				break;
			}

			sch.erase(attr);
		}

		attr = sch.find("default");
		if (attr != sch.end()) {
			set_default_value(attr.value());
			sch.erase(attr);
		}

		for (auto &key : known_keywords)
			sch.erase(key);

		// with nlohmann::json float instance (but number in schema-definition) can be seen as unsigned or integer -
		// reuse the number-validator for integer values as well, if they have not been specified explicitly
		if (type_[static_cast<uint8_t>(json::value_t::number_float)] && !type_[static_cast<uint8_t>(json::value_t::number_integer)])
			type_[static_cast<uint8_t>(json::value_t::number_integer)] = type_[static_cast<uint8_t>(json::value_t::number_float)];

		// #54: JSON-schema does not differentiate between unsigned and signed integer - nlohmann::json does
		// we stick with JSON-schema: use the integer-validator if instance-value is unsigned
		type_[static_cast<uint8_t>(json::value_t::number_unsigned)] = type_[static_cast<uint8_t>(json::value_t::number_integer)];

		// special for binary types
		if (type_[static_cast<uint8_t>(json::value_t::string)]) {
			type_[static_cast<uint8_t>(json::value_t::binary)] = type_[static_cast<uint8_t>(json::value_t::string)];
		}

		attr = sch.find("enum");
		if (attr != sch.end()) {
			enum_ = {true, attr.value()};
			sch.erase(attr);
		}

		attr = sch.find("const");
		if (attr != sch.end()) {
			const_ = {true, attr.value()};
			sch.erase(attr);
		}

		attr = sch.find("not");
		if (attr != sch.end()) {
			logic_.push_back(std::make_shared<logical_not>(attr.value(), root, uris));
			sch.erase(attr);
		}

		attr = sch.find("allOf");
		if (attr != sch.end()) {
			logic_.push_back(std::make_shared<logical_combination<allOf>>(attr.value(), root, uris));
			sch.erase(attr);
		}

		attr = sch.find("anyOf");
		if (attr != sch.end()) {
			logic_.push_back(std::make_shared<logical_combination<anyOf>>(attr.value(), root, uris));
			sch.erase(attr);
		}

		attr = sch.find("oneOf");
		if (attr != sch.end()) {
			logic_.push_back(std::make_shared<logical_combination<oneOf>>(attr.value(), root, uris));
			sch.erase(attr);
		}

		attr = sch.find("if");
		if (attr != sch.end()) {
			auto attr_then = sch.find("then");
			auto attr_else = sch.find("else");

			if (attr_then != sch.end() || attr_else != sch.end()) {
				if_ = schema::make(attr.value(), root, {"if"}, uris);

				if (attr_then != sch.end()) {
					then_ = schema::make(attr_then.value(), root, {"then"}, uris);
					sch.erase(attr_then);
				}

				if (attr_else != sch.end()) {
					else_ = schema::make(attr_else.value(), root, {"else"}, uris);
					sch.erase(attr_else);
				}
			}
			sch.erase(attr);
		}
	}
};

class string : public schema
{
	std::pair<bool, size_t> maxLength_{false, 0};
	std::pair<bool, size_t> minLength_{false, 0};

#ifndef NO_STD_REGEX
	std::pair<bool, REGEX_NAMESPACE::regex> pattern_{false, REGEX_NAMESPACE::regex()};
	std::string patternString_;
#endif

	std::pair<bool, std::string> format_;
	std::tuple<bool, std::string, std::string> content_{false, "", ""};

	std::size_t utf8_length(const std::string &s) const
	{
		size_t len = 0;
		for (auto c : s)
			if ((c & 0xc0) != 0x80)
				len++;
		return len;
	}

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &, error_handler &e) const override
	{
		if (minLength_.first) {
			if (utf8_length(instance.get<std::string>()) < minLength_.second) {
				std::ostringstream s;
				s << "instance is too short as per minLength:" << minLength_.second;
				e.error(ptr, instance, s.str());
			}
		}

		if (maxLength_.first) {
			if (utf8_length(instance.get<std::string>()) > maxLength_.second) {
				std::ostringstream s;
				s << "instance is too long as per maxLength: " << maxLength_.second;
				e.error(ptr, instance, s.str());
			}
		}

		if (std::get<0>(content_)) {
			if (root_->content_check() == nullptr)
				e.error(ptr, instance, std::string("a content checker was not provided but a contentEncoding or contentMediaType for this string have been present: '") + std::get<1>(content_) + "' '" + std::get<2>(content_) + "'");
			else {
				try {
					root_->content_check()(std::get<1>(content_), std::get<2>(content_), instance);
				} catch (const std::exception &ex) {
					e.error(ptr, instance, std::string("content-checking failed: ") + ex.what());
				}
			}
		} else if (instance.type() == json::value_t::binary) {
			e.error(ptr, instance, "expected string, but get binary data");
		}

		if (instance.type() != json::value_t::string) {
			return; // next checks only for strings
		}

#ifndef NO_STD_REGEX
		if (pattern_.first &&
		    !REGEX_NAMESPACE::regex_search(instance.get<std::string>(), pattern_.second))
			e.error(ptr, instance, "instance does not match regex pattern: " + patternString_);
#endif

		if (format_.first) {
			if (root_->format_check() == nullptr)
				e.error(ptr, instance, std::string("a format checker was not provided but a format keyword for this string is present: ") + format_.second);
			else {
				try {
					root_->format_check()(format_.second, instance.get<std::string>());
				} catch (const std::exception &ex) {
					e.error(ptr, instance, std::string("format-checking failed: ") + ex.what());
				}
			}
		}
	}

public:
	string(json &sch, root_schema *root)
	    : schema(root)
	{
		auto attr = sch.find("maxLength");
		if (attr != sch.end()) {
			maxLength_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("minLength");
		if (attr != sch.end()) {
			minLength_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("contentEncoding");
		if (attr != sch.end()) {
			std::get<0>(content_) = true;
			std::get<1>(content_) = attr.value().get<std::string>();

			// special case for nlohmann::json-binary-types
			//
			// https://github.com/pboettch/json-schema-validator/pull/114
			//
			// We cannot use explicitly in a schema: {"type": "binary"} or
			// "type": ["binary", "number"] we have to be implicit. For a
			// schema where "contentEncoding" is set to "binary", an instance
			// of type json::value_t::binary is accepted. If a
			// contentEncoding-callback has to be provided and is called
			// accordingly. For encoding=binary, no other type validations are done

			sch.erase(attr);
		}

		attr = sch.find("contentMediaType");
		if (attr != sch.end()) {
			std::get<0>(content_) = true;
			std::get<2>(content_) = attr.value().get<std::string>();

			sch.erase(attr);
		}

		if (std::get<0>(content_) == true && root_->content_check() == nullptr) {
			throw std::invalid_argument{"schema contains contentEncoding/contentMediaType but content checker was not set"};
		}

#ifndef NO_STD_REGEX
		attr = sch.find("pattern");
		if (attr != sch.end()) {
			patternString_ = attr.value().get<std::string>();
			pattern_ = {true, REGEX_NAMESPACE::regex(attr.value().get<std::string>(),
			                                         REGEX_NAMESPACE::regex::ECMAScript)};
			sch.erase(attr);
		}
#endif

		attr = sch.find("format");
		if (attr != sch.end()) {
			if (root_->format_check() == nullptr)
				throw std::invalid_argument{"a format checker was not provided but a format keyword for this string is present: " + format_.second};

			format_ = {true, attr.value().get<std::string>()};
			sch.erase(attr);
		}
	}
};

template <typename T>
class numeric : public schema
{
	std::pair<bool, T> maximum_{false, 0};
	std::pair<bool, T> minimum_{false, 0};

	bool exclusiveMaximum_ = false;
	bool exclusiveMinimum_ = false;

	std::pair<bool, json::number_float_t> multipleOf_{false, 0};

	// multipleOf - if the remainder of the division is 0 -> OK
	bool violates_multiple_of(T x) const
	{
		double res = std::remainder(x, multipleOf_.second);
		double multiple = std::fabs(x / multipleOf_.second);
		if (multiple > 1) {
			res = res / multiple;
		}
		double eps = std::nextafter(x, 0) - static_cast<double>(x);

		return std::fabs(res) > std::fabs(eps);
	}

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &, error_handler &e) const override
	{
		T value = instance; // conversion of json to value_type

		std::ostringstream oss;

		if (multipleOf_.first && value != 0) // zero is multiple of everything
			if (violates_multiple_of(value))
				oss << "instance is not a multiple of " << json(multipleOf_.second);

		if (maximum_.first) {
			if (exclusiveMaximum_ && value >= maximum_.second)
				oss << "instance exceeds or equals maximum of " << json(maximum_.second);
			else if (value > maximum_.second)
				oss << "instance exceeds maximum of " << json(maximum_.second);
		}

		if (minimum_.first) {
			if (exclusiveMinimum_ && value <= minimum_.second)
				oss << "instance is below or equals minimum of " << json(minimum_.second);
			else if (value < minimum_.second)
				oss << "instance is below minimum of " << json(minimum_.second);
		}

		oss.seekp(0, std::ios::end);
		auto size = oss.tellp();
		if (size != 0) {
			oss.seekp(0, std::ios::beg);
			e.error(ptr, instance, oss.str());
		}
	}

public:
	numeric(const json &sch, root_schema *root, std::set<std::string> &kw)
	    : schema(root)
	{
		auto attr = sch.find("maximum");
		if (attr != sch.end()) {
			maximum_ = {true, attr.value().get<T>()};
			kw.insert("maximum");
		}

		attr = sch.find("minimum");
		if (attr != sch.end()) {
			minimum_ = {true, attr.value().get<T>()};
			kw.insert("minimum");
		}

		attr = sch.find("exclusiveMaximum");
		if (attr != sch.end()) {
			exclusiveMaximum_ = true;
			maximum_ = {true, attr.value().get<T>()};
			kw.insert("exclusiveMaximum");
		}

		attr = sch.find("exclusiveMinimum");
		if (attr != sch.end()) {
			exclusiveMinimum_ = true;
			minimum_ = {true, attr.value().get<T>()};
			kw.insert("exclusiveMinimum");
		}

		attr = sch.find("multipleOf");
		if (attr != sch.end()) {
			multipleOf_ = {true, attr.value().get<json::number_float_t>()};
			kw.insert("multipleOf");
		}
	}
};

class null : public schema
{
	void validate(const json::json_pointer &ptr, const json &instance, json_patch &, error_handler &e) const override
	{
		if (!instance.is_null())
			e.error(ptr, instance, "expected to be null");
	}

public:
	null(json &, root_schema *root)
	    : schema(root) {}
};

class boolean_type : public schema
{
	void validate(const json::json_pointer &, const json &, json_patch &, error_handler &) const override {}

public:
	boolean_type(json &, root_schema *root)
	    : schema(root) {}
};

class boolean : public schema
{
	bool true_;
	void validate(const json::json_pointer &ptr, const json &instance, json_patch &, error_handler &e) const override
	{
		if (!true_) { // false schema
			// empty array
			// switch (instance.type()) {
			// case json::value_t::array:
			//	if (instance.size() != 0) // valid false-schema
			//		e.error(ptr, instance, "false-schema required empty array");
			//	return;
			//}

			e.error(ptr, instance, "instance invalid as per false-schema");
		}
	}

public:
	boolean(json &sch, root_schema *root)
	    : schema(root), true_(sch) {}
};

class required : public schema
{
	const std::vector<std::string> required_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &, error_handler &e) const override final
	{
		for (auto &r : required_)
			if (instance.find(r) == instance.end())
				e.error(ptr, instance, "required property '" + r + "' not found in object as a dependency");
	}

public:
	required(const std::vector<std::string> &r, root_schema *root)
	    : schema(root), required_(r) {}
};

class object : public schema
{
	std::pair<bool, size_t> maxProperties_{false, 0};
	std::pair<bool, size_t> minProperties_{false, 0};
	std::vector<std::string> required_;

	std::map<std::string, std::shared_ptr<schema>> properties_;
#ifndef NO_STD_REGEX
	std::vector<std::pair<REGEX_NAMESPACE::regex, std::shared_ptr<schema>>> patternProperties_;
#endif
	std::shared_ptr<schema> additionalProperties_;

	std::map<std::string, std::shared_ptr<schema>> dependencies_;

	std::shared_ptr<schema> propertyNames_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const override
	{
		if (maxProperties_.first && instance.size() > maxProperties_.second)
			e.error(ptr, instance, "too many properties");

		if (minProperties_.first && instance.size() < minProperties_.second)
			e.error(ptr, instance, "too few properties");

		for (auto &r : required_)
			if (instance.find(r) == instance.end())
				e.error(ptr, instance, "required property '" + r + "' not found in object");

		// for each property in instance
		for (auto &p : instance.items()) {
			if (propertyNames_)
				propertyNames_->validate(ptr, p.key(), patch, e);

			bool a_prop_or_pattern_matched = false;
			auto schema_p = properties_.find(p.key());
			// check if it is in "properties"
			if (schema_p != properties_.end()) {
				a_prop_or_pattern_matched = true;
				schema_p->second->validate(ptr / p.key(), p.value(), patch, e);
			}

#ifndef NO_STD_REGEX
			// check all matching patternProperties
			for (auto &schema_pp : patternProperties_)
				if (REGEX_NAMESPACE::regex_search(p.key(), schema_pp.first)) {
					a_prop_or_pattern_matched = true;
					schema_pp.second->validate(ptr / p.key(), p.value(), patch, e);
				}
#endif

			// check additionalProperties as a last resort
			if (!a_prop_or_pattern_matched && additionalProperties_) {
				first_error_handler additional_prop_err;
				additionalProperties_->validate(ptr / p.key(), p.value(), patch, additional_prop_err);
				if (additional_prop_err)
					e.error(ptr, instance, "validation failed for additional property '" + p.key() + "': " + additional_prop_err.message_);
			}
		}

		// reverse search
		for (auto const &prop : properties_) {
			const auto finding = instance.find(prop.first);
			if (instance.end() == finding) { // if the prop is not in the instance
				const auto &default_value = prop.second->default_value(ptr, instance, e);
				if (!default_value.is_null()) { // if default value is available
					patch.add((ptr / prop.first), default_value);
				}
			}
		}

		for (auto &dep : dependencies_) {
			auto prop = instance.find(dep.first);
			if (prop != instance.end())                                    // if dependency-property is present in instance
				dep.second->validate(ptr / dep.first, instance, patch, e); // validate
		}
	}

public:
	object(json &sch,
	       root_schema *root,
	       const std::vector<nlohmann::json_uri> &uris)
	    : schema(root)
	{
		auto attr = sch.find("maxProperties");
		if (attr != sch.end()) {
			maxProperties_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("minProperties");
		if (attr != sch.end()) {
			minProperties_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("required");
		if (attr != sch.end()) {
			required_ = attr.value().get<std::vector<std::string>>();
			sch.erase(attr);
		}

		attr = sch.find("properties");
		if (attr != sch.end()) {
			for (auto prop : attr.value().items())
				properties_.insert(
				    std::make_pair(
				        prop.key(),
				        schema::make(prop.value(), root, {"properties", prop.key()}, uris)));
			sch.erase(attr);
		}

#ifndef NO_STD_REGEX
		attr = sch.find("patternProperties");
		if (attr != sch.end()) {
			for (auto prop : attr.value().items())
				patternProperties_.push_back(
				    std::make_pair(
				        REGEX_NAMESPACE::regex(prop.key(), REGEX_NAMESPACE::regex::ECMAScript),
				        schema::make(prop.value(), root, {prop.key()}, uris)));
			sch.erase(attr);
		}
#endif

		attr = sch.find("additionalProperties");
		if (attr != sch.end()) {
			additionalProperties_ = schema::make(attr.value(), root, {"additionalProperties"}, uris);
			sch.erase(attr);
		}

		attr = sch.find("dependencies");
		if (attr != sch.end()) {
			for (auto &dep : attr.value().items())
				switch (dep.value().type()) {
				case json::value_t::array:
					dependencies_.emplace(dep.key(),
					                      std::make_shared<required>(
					                          dep.value().get<std::vector<std::string>>(), root));
					break;

				default:
					dependencies_.emplace(dep.key(),
					                      schema::make(dep.value(), root, {"dependencies", dep.key()}, uris));
					break;
				}
			sch.erase(attr);
		}

		attr = sch.find("propertyNames");
		if (attr != sch.end()) {
			propertyNames_ = schema::make(attr.value(), root, {"propertyNames"}, uris);
			sch.erase(attr);
		}

		attr = sch.find("default");
		if (attr != sch.end()) {
			set_default_value(*attr);
		}
	}
};

class array : public schema
{
	std::pair<bool, size_t> maxItems_{false, 0};
	std::pair<bool, size_t> minItems_{false, 0};
	bool uniqueItems_ = false;

	std::shared_ptr<schema> items_schema_;

	std::vector<std::shared_ptr<schema>> items_;
	std::shared_ptr<schema> additionalItems_;

	std::shared_ptr<schema> contains_;

	void validate(const json::json_pointer &ptr, const json &instance, json_patch &patch, error_handler &e) const override
	{
		if (maxItems_.first && instance.size() > maxItems_.second)
			e.error(ptr, instance, "array has too many items");

		if (minItems_.first && instance.size() < minItems_.second)
			e.error(ptr, instance, "array has too few items");

		if (uniqueItems_) {
			for (auto it = instance.cbegin(); it != instance.cend(); ++it) {
				auto v = std::find(it + 1, instance.end(), *it);
				if (v != instance.end())
					e.error(ptr, instance, "items have to be unique for this array");
			}
		}

		size_t index = 0;
		if (items_schema_)
			for (auto &i : instance) {
				items_schema_->validate(ptr / index, i, patch, e);
				index++;
			}
		else {
			auto item = items_.cbegin();
			for (auto &i : instance) {
				std::shared_ptr<schema> item_validator;
				if (item == items_.cend())
					item_validator = additionalItems_;
				else {
					item_validator = *item;
					item++;
				}

				if (!item_validator)
					break;

				item_validator->validate(ptr / index, i, patch, e);
			}
		}

		if (contains_) {
			bool contained = false;
			for (auto &item : instance) {
				first_error_handler local_e;
				contains_->validate(ptr, item, patch, local_e);
				if (!local_e) {
					contained = true;
					break;
				}
			}
			if (!contained)
				e.error(ptr, instance, "array does not contain required element as per 'contains'");
		}
	}

public:
	array(json &sch, root_schema *root, const std::vector<nlohmann::json_uri> &uris)
	    : schema(root)
	{
		auto attr = sch.find("maxItems");
		if (attr != sch.end()) {
			maxItems_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("minItems");
		if (attr != sch.end()) {
			minItems_ = {true, attr.value().get<size_t>()};
			sch.erase(attr);
		}

		attr = sch.find("uniqueItems");
		if (attr != sch.end()) {
			uniqueItems_ = attr.value().get<bool>();
			sch.erase(attr);
		}

		attr = sch.find("items");
		if (attr != sch.end()) {

			if (attr.value().type() == json::value_t::array) {
				size_t c = 0;
				for (auto &subsch : attr.value())
					items_.push_back(schema::make(subsch, root, {"items", std::to_string(c++)}, uris));

				auto attr_add = sch.find("additionalItems");
				if (attr_add != sch.end()) {
					additionalItems_ = schema::make(attr_add.value(), root, {"additionalItems"}, uris);
					sch.erase(attr_add);
				}

			} else if (attr.value().type() == json::value_t::object ||
			           attr.value().type() == json::value_t::boolean)
				items_schema_ = schema::make(attr.value(), root, {"items"}, uris);

			sch.erase(attr);
		}

		attr = sch.find("contains");
		if (attr != sch.end()) {
			contains_ = schema::make(attr.value(), root, {"contains"}, uris);
			sch.erase(attr);
		}
	}
};

std::shared_ptr<schema> type_schema::make(json &schema,
                                          json::value_t type,
                                          root_schema *root,
                                          const std::vector<nlohmann::json_uri> &uris,
                                          std::set<std::string> &kw)
{
	switch (type) {
	case json::value_t::null:
		return std::make_shared<null>(schema, root);

	case json::value_t::number_unsigned:
	case json::value_t::number_integer:
		return std::make_shared<numeric<json::number_integer_t>>(schema, root, kw);
	case json::value_t::number_float:
		return std::make_shared<numeric<json::number_float_t>>(schema, root, kw);
	case json::value_t::string:
		return std::make_shared<string>(schema, root);
	case json::value_t::boolean:
		return std::make_shared<boolean_type>(schema, root);
	case json::value_t::object:
		return std::make_shared<object>(schema, root, uris);
	case json::value_t::array:
		return std::make_shared<array>(schema, root, uris);

	case json::value_t::discarded: // not a real type - silence please
		break;

	case json::value_t::binary:
		break;
	}
	return nullptr;
}
} // namespace

namespace
{

std::shared_ptr<schema> schema::make(json &schema,
                                     root_schema *root,
                                     const std::vector<std::string> &keys,
                                     std::vector<nlohmann::json_uri> uris)
{
	// remove URIs which contain plain name identifiers, as sub-schemas cannot be referenced
	for (auto uri = uris.begin(); uri != uris.end();)
		if (uri->identifier() != "")
			uri = uris.erase(uri);
		else
			uri++;

	// append to all URIs the keys for this sub-schema
	for (auto &key : keys)
		for (auto &uri : uris)
			uri = uri.append(key);

	std::shared_ptr<::schema> sch;

	// boolean schema
	if (schema.type() == json::value_t::boolean)
		sch = std::make_shared<boolean>(schema, root);
	else if (schema.type() == json::value_t::object) {

		auto attr = schema.find("$id"); // if $id is present, this schema can be referenced by this ID
		                                // as an additional URI
		if (attr != schema.end()) {
			if (std::find(uris.begin(),
			              uris.end(),
			              attr.value().get<std::string>()) == uris.end())
				uris.push_back(uris.back().derive(attr.value().get<std::string>())); // so add it to the list if it is not there already
			schema.erase(attr);
		}

		auto findDefinitions = [&](const std::string &defs) -> bool {
			attr = schema.find(defs);
			if (attr != schema.end()) {
				for (auto &def : attr.value().items())
					schema::make(def.value(), root, {defs, def.key()}, uris);
				schema.erase(attr);
				return true;
			}
			return false;
		};
		if (!findDefinitions("$defs")) {
			findDefinitions("definitions");
		}

		attr = schema.find("$ref");
		if (attr != schema.end()) { // this schema is a reference
			// the last one on the uri-stack is the last id seen before coming here,
			// so this is the origial URI for this reference, the $ref-value has thus be resolved from it
			auto id = uris.back().derive(attr.value().get<std::string>());
			sch = root->get_or_create_ref(id);

			schema.erase(attr);

			// special case where we break draft-7 and allow overriding of properties when a $ref is used
			attr = schema.find("default");
			if (attr != schema.end()) {
				// copy the referenced schema depending on the underlying type and modify the default value
				if (auto new_sch = sch->make_for_default_(sch, root, uris, attr.value())) {
					sch = new_sch;
				}
				schema.erase(attr);
			}
		} else {
			sch = std::make_shared<type_schema>(schema, root, uris);
		}

		schema.erase("$schema");
		schema.erase("title");
		schema.erase("description");
	} else {
		throw std::invalid_argument("invalid JSON-type for a schema for " + uris[0].to_string() + ", expected: boolean or object");
	}

	for (auto &uri : uris) { // for all URIs this schema is referenced by
		root->insert(uri, sch);

		if (schema.type() == json::value_t::object)
			for (auto &u : schema.items())
				root->insert_unknown_keyword(uri, u.key(), u.value()); // insert unknown keywords for later reference
	}
	return sch;
}

class throwing_error_handler : public error_handler
{
	void error(const json::json_pointer &ptr, const json &instance, const std::string &message) override
	{
		throw std::invalid_argument(std::string("At ") + ptr.to_string() + " of " + instance.dump() + " - " + message + "\n");
	}
};

} // namespace

namespace nlohmann
{
namespace json_schema
{

json_validator::json_validator(schema_loader loader,
                               format_checker format,
                               content_checker content)
    : root_(std::unique_ptr<root_schema>(new root_schema(std::move(loader),
                                                         std::move(format),
                                                         std::move(content))))
{
}

json_validator::json_validator(const json &schema,
                               schema_loader loader,
                               format_checker format,
                               content_checker content)
    : json_validator(std::move(loader),
                     std::move(format),
                     std::move(content))
{
	set_root_schema(schema);
}

json_validator::json_validator(json &&schema,
                               schema_loader loader,
                               format_checker format,
                               content_checker content)

    : json_validator(std::move(loader),
                     std::move(format),
                     std::move(content))
{
	set_root_schema(std::move(schema));
}

// move constructor, destructor and move assignment operator can be defaulted here
// where root_schema is a complete type
json_validator::json_validator(json_validator &&) = default;
json_validator::~json_validator() = default;
json_validator &json_validator::operator=(json_validator &&) = default;

void json_validator::set_root_schema(const json &schema)
{
	root_->set_root_schema(schema);
}

void json_validator::set_root_schema(json &&schema)
{
	root_->set_root_schema(std::move(schema));
}

json json_validator::validate(const json &instance) const
{
	throwing_error_handler err;
	return validate(instance, err);
}

json json_validator::validate(const json &instance, error_handler &err, const json_uri &initial_uri) const
{
	json::json_pointer ptr;
	json_patch patch;
	root_->validate(ptr, instance, patch, err, initial_uri);
	return patch;
}

} // namespace json_schema
} // namespace nlohmann
#ifndef SMTP_ADDRESS_PARSER_HPP_INCLUDED
#define SMTP_ADDRESS_PARSER_HPP_INCLUDED

/*

Snarfed from <https://github.com/gene-hightower/smtp-address-validator>

<http://opensource.org/licenses/MIT>:

Copyright (c) 2021 Gene Hightower

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

bool is_address(const char *p, const char *pe);

#endif // SMTP_ADDRESS_PARSER_HPP_INCLUDED
/*

Snarfed from <https://github.com/gene-hightower/smtp-address-validator>

<http://opensource.org/licenses/MIT>:

Copyright (c) 2021 Gene Hightower

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

// #include "smtp-address-validator.hpp"
#ifndef SMTP_ADDRESS_PARSER_HPP_INCLUDED
#define SMTP_ADDRESS_PARSER_HPP_INCLUDED

/*

Snarfed from <https://github.com/gene-hightower/smtp-address-validator>

<http://opensource.org/licenses/MIT>:

Copyright (c) 2021 Gene Hightower

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

bool is_address(const char *p, const char *pe);

#endif // SMTP_ADDRESS_PARSER_HPP_INCLUDED


static const signed char _address_actions[] = {
    0, 1, 0, 1, 1, 0};

static const short _address_key_offsets[] = {
    0, 0, 24, 26, 50, 52, 54, 56,
    58, 60, 62, 86, 103, 105, 107, 109,
    111, 113, 115, 117, 134, 150, 161, 168,
    176, 180, 181, 190, 195, 196, 201, 202,
    207, 210, 213, 219, 222, 225, 228, 234,
    237, 240, 243, 249, 252, 261, 270, 282,
    293, 302, 311, 320, 328, 345, 353, 360,
    367, 368, 375, 382, 389, 396, 397, 404,
    411, 418, 425, 426, 433, 440, 447, 454,
    455, 462, 469, 476, 483, 484, 491, 498,
    505, 512, 513, 523, 531, 538, 545, 546,
    552, 559, 566, 573, 581, 589, 597, 608,
    618, 626, 634, 641, 649, 657, 665, 667,
    673, 681, 689, 697, 699, 705, 713, 721,
    729, 731, 737, 745, 753, 761, 763, 769,
    777, 785, 793, 795, 802, 812, 821, 829,
    837, 839, 848, 857, 865, 873, 875, 884,
    893, 901, 909, 911, 920, 929, 937, 945,
    947, 956, 965, 974, 983, 992, 1004, 1015,
    1024, 1033, 1042, 1051, 1060, 1072, 1083, 1092,
    1101, 1109, 1118, 1127, 1136, 1148, 1159, 1168,
    1177, 1185, 1194, 1203, 1212, 1224, 1235, 1244,
    1253, 1261, 1270, 1279, 1288, 1300, 1311, 1320,
    1329, 1337, 1339, 1353, 1355, 1357, 1359, 1361,
    1363, 1365, 1367, 1368, 1370, 1388, 0};

static const signed char _address_trans_keys[] = {
    -32, -19, -16, -12, 34, 45, 61, 63,
    -62, -33, -31, -17, -15, -13, 33, 39,
    42, 43, 47, 57, 65, 90, 94, 126,
    -128, -65, -32, -19, -16, -12, 33, 46,
    61, 64, -62, -33, -31, -17, -15, -13,
    35, 39, 42, 43, 45, 57, 63, 90,
    94, 126, -96, -65, -128, -65, -128, -97,
    -112, -65, -128, -65, -128, -113, -32, -19,
    -16, -12, 33, 45, 61, 63, -62, -33,
    -31, -17, -15, -13, 35, 39, 42, 43,
    47, 57, 65, 90, 94, 126, -32, -19,
    -16, -12, 91, -62, -33, -31, -17, -15,
    -13, 48, 57, 65, 90, 97, 122, -128,
    -65, -96, -65, -128, -65, -128, -97, -112,
    -65, -128, -65, -128, -113, -32, -19, -16,
    -12, 45, -62, -33, -31, -17, -15, -13,
    48, 57, 65, 90, 97, 122, -32, -19,
    -16, -12, -62, -33, -31, -17, -15, -13,
    48, 57, 65, 90, 97, 122, 45, 48,
    49, 50, 73, 51, 57, 65, 90, 97,
    122, 45, 48, 57, 65, 90, 97, 122,
    45, 58, 48, 57, 65, 90, 97, 122,
    33, 90, 94, 126, 93, 45, 46, 58,
    48, 57, 65, 90, 97, 122, 48, 49,
    50, 51, 57, 46, 48, 49, 50, 51,
    57, 46, 48, 49, 50, 51, 57, 93,
    48, 57, 93, 48, 57, 53, 93, 48,
    52, 54, 57, 93, 48, 53, 46, 48,
    57, 46, 48, 57, 46, 53, 48, 52,
    54, 57, 46, 48, 53, 46, 48, 57,
    46, 48, 57, 46, 53, 48, 52, 54,
    57, 46, 48, 53, 45, 46, 58, 48,
    57, 65, 90, 97, 122, 45, 46, 58,
    48, 57, 65, 90, 97, 122, 45, 46,
    53, 58, 48, 52, 54, 57, 65, 90,
    97, 122, 45, 46, 58, 48, 53, 54,
    57, 65, 90, 97, 122, 45, 58, 80,
    48, 57, 65, 90, 97, 122, 45, 58,
    118, 48, 57, 65, 90, 97, 122, 45,
    54, 58, 48, 57, 65, 90, 97, 122,
    45, 58, 48, 57, 65, 90, 97, 122,
    58, 33, 47, 48, 57, 59, 64, 65,
    70, 71, 90, 94, 96, 97, 102, 103,
    126, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 48, 57, 65, 70, 97, 102,
    58, 48, 57, 65, 70, 97, 102, 58,
    58, 48, 57, 65, 70, 97, 102, 58,
    48, 57, 65, 70, 97, 102, 58, 48,
    57, 65, 70, 97, 102, 58, 48, 57,
    65, 70, 97, 102, 58, 58, 48, 57,
    65, 70, 97, 102, 58, 48, 57, 65,
    70, 97, 102, 58, 48, 57, 65, 70,
    97, 102, 58, 48, 57, 65, 70, 97,
    102, 58, 58, 48, 57, 65, 70, 97,
    102, 58, 48, 57, 65, 70, 97, 102,
    58, 48, 57, 65, 70, 97, 102, 58,
    48, 57, 65, 70, 97, 102, 58, 58,
    48, 57, 65, 70, 97, 102, 58, 48,
    57, 65, 70, 97, 102, 58, 48, 57,
    65, 70, 97, 102, 58, 48, 57, 65,
    70, 97, 102, 58, 58, 48, 57, 65,
    70, 97, 102, 58, 48, 57, 65, 70,
    97, 102, 58, 48, 57, 65, 70, 97,
    102, 58, 48, 57, 65, 70, 97, 102,
    58, 48, 49, 50, 58, 51, 57, 65,
    70, 97, 102, 46, 58, 48, 57, 65,
    70, 97, 102, 58, 48, 57, 65, 70,
    97, 102, 58, 48, 57, 65, 70, 97,
    102, 58, 48, 57, 65, 70, 97, 102,
    93, 48, 57, 65, 70, 97, 102, 93,
    48, 57, 65, 70, 97, 102, 93, 48,
    57, 65, 70, 97, 102, 46, 58, 48,
    57, 65, 70, 97, 102, 46, 58, 48,
    57, 65, 70, 97, 102, 46, 58, 48,
    57, 65, 70, 97, 102, 46, 53, 58,
    48, 52, 54, 57, 65, 70, 97, 102,
    46, 58, 48, 53, 54, 57, 65, 70,
    97, 102, 46, 58, 48, 57, 65, 70,
    97, 102, 46, 58, 48, 57, 65, 70,
    97, 102, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 58, 48, 57, 65, 70,
    97, 102, 48, 49, 50, 93, 51, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    49, 50, 51, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 49, 50, 51, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    49, 50, 51, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 49, 50, 51, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    57, 65, 70, 97, 102, 46, 58, 93,
    48, 57, 65, 70, 97, 102, 46, 58,
    93, 48, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 53, 58, 93, 48, 52, 54, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    53, 54, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 57, 65, 70, 97,
    102, 46, 58, 93, 48, 57, 65, 70,
    97, 102, 46, 58, 93, 48, 57, 65,
    70, 97, 102, 46, 58, 93, 48, 57,
    65, 70, 97, 102, 46, 53, 58, 93,
    48, 52, 54, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 53, 54, 57, 65,
    70, 97, 102, 46, 58, 93, 48, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 46, 58, 93,
    48, 57, 65, 70, 97, 102, 46, 58,
    93, 48, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 53, 58, 93, 48, 52, 54, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    53, 54, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 46, 58, 93, 48, 57, 65, 70,
    97, 102, 46, 58, 93, 48, 57, 65,
    70, 97, 102, 46, 58, 93, 48, 57,
    65, 70, 97, 102, 46, 53, 58, 93,
    48, 52, 54, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 53, 54, 57, 65,
    70, 97, 102, 46, 58, 93, 48, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    57, 65, 70, 97, 102, 58, 93, 48,
    57, 65, 70, 97, 102, 46, 58, 93,
    48, 57, 65, 70, 97, 102, 46, 58,
    93, 48, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 53, 58, 93, 48, 52, 54, 57,
    65, 70, 97, 102, 46, 58, 93, 48,
    53, 54, 57, 65, 70, 97, 102, 46,
    58, 93, 48, 57, 65, 70, 97, 102,
    46, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, 48, 57, 65, 70, 97,
    102, 58, 93, -32, -19, -16, -12, 34,
    92, -62, -33, -31, -17, -15, -13, 32,
    126, -128, -65, -96, -65, -128, -65, -128,
    -97, -112, -65, -128, -65, -128, -113, 64,
    32, 126, -32, -19, -16, -12, 45, 46,
    -62, -33, -31, -17, -15, -13, 48, 57,
    65, 90, 97, 122, 0};

static const signed char _address_single_lengths[] = {
    0, 8, 0, 8, 0, 0, 0, 0,
    0, 0, 8, 5, 0, 0, 0, 0,
    0, 0, 0, 5, 4, 5, 1, 2,
    0, 1, 3, 3, 1, 3, 1, 3,
    1, 1, 2, 1, 1, 1, 2, 1,
    1, 1, 2, 1, 3, 3, 4, 3,
    3, 3, 3, 2, 1, 2, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 4, 2, 1, 1, 1, 0,
    1, 1, 1, 2, 2, 2, 3, 2,
    2, 2, 1, 2, 2, 2, 2, 0,
    2, 2, 2, 2, 0, 2, 2, 2,
    2, 0, 2, 2, 2, 2, 0, 2,
    2, 2, 2, 1, 4, 3, 2, 2,
    2, 3, 3, 2, 2, 2, 3, 3,
    2, 2, 2, 3, 3, 2, 2, 2,
    3, 3, 3, 3, 3, 4, 3, 3,
    3, 3, 3, 3, 4, 3, 3, 3,
    2, 3, 3, 3, 4, 3, 3, 3,
    2, 3, 3, 3, 4, 3, 3, 3,
    2, 3, 3, 3, 4, 3, 3, 3,
    2, 2, 6, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 6, 0, 0};

static const signed char _address_range_lengths[] = {
    0, 8, 1, 8, 1, 1, 1, 1,
    1, 1, 8, 6, 1, 1, 1, 1,
    1, 1, 1, 6, 6, 3, 3, 3,
    2, 0, 3, 1, 0, 1, 0, 1,
    1, 1, 2, 1, 1, 1, 2, 1,
    1, 1, 2, 1, 3, 3, 4, 4,
    3, 3, 3, 3, 8, 3, 3, 3,
    0, 3, 3, 3, 3, 0, 3, 3,
    3, 3, 0, 3, 3, 3, 3, 0,
    3, 3, 3, 3, 0, 3, 3, 3,
    3, 0, 3, 3, 3, 3, 0, 3,
    3, 3, 3, 3, 3, 3, 4, 4,
    3, 3, 3, 3, 3, 3, 0, 3,
    3, 3, 3, 0, 3, 3, 3, 3,
    0, 3, 3, 3, 3, 0, 3, 3,
    3, 3, 0, 3, 3, 3, 3, 3,
    0, 3, 3, 3, 3, 0, 3, 3,
    3, 3, 0, 3, 3, 3, 3, 0,
    3, 3, 3, 3, 3, 4, 4, 3,
    3, 3, 3, 3, 4, 4, 3, 3,
    3, 3, 3, 3, 4, 4, 3, 3,
    3, 3, 3, 3, 4, 4, 3, 3,
    3, 3, 3, 3, 4, 4, 3, 3,
    3, 0, 4, 1, 1, 1, 1, 1,
    1, 1, 0, 1, 6, 0, 0};

static const short _address_index_offsets[] = {
    0, 0, 17, 19, 36, 38, 40, 42,
    44, 46, 48, 65, 77, 79, 81, 83,
    85, 87, 89, 91, 103, 114, 123, 128,
    134, 137, 139, 146, 151, 153, 158, 160,
    165, 168, 171, 176, 179, 182, 185, 190,
    193, 196, 199, 204, 207, 214, 221, 230,
    238, 245, 252, 259, 265, 275, 281, 286,
    291, 293, 298, 303, 308, 313, 315, 320,
    325, 330, 335, 337, 342, 347, 352, 357,
    359, 364, 369, 374, 379, 381, 386, 391,
    396, 401, 403, 411, 417, 422, 427, 429,
    433, 438, 443, 448, 454, 460, 466, 474,
    481, 487, 493, 498, 504, 510, 516, 519,
    523, 529, 535, 541, 544, 548, 554, 560,
    566, 569, 573, 579, 585, 591, 594, 598,
    604, 610, 616, 619, 624, 632, 639, 645,
    651, 654, 661, 668, 674, 680, 683, 690,
    697, 703, 709, 712, 719, 726, 732, 738,
    741, 748, 755, 762, 769, 776, 785, 793,
    800, 807, 814, 821, 828, 837, 845, 852,
    859, 865, 872, 879, 886, 895, 903, 910,
    917, 923, 930, 937, 944, 953, 961, 968,
    975, 981, 988, 995, 1002, 1011, 1019, 1026,
    1033, 1039, 1042, 1053, 1055, 1057, 1059, 1061,
    1063, 1065, 1067, 1069, 1071, 1084, 0};

static const short _address_cond_targs[] = {
    4, 6, 7, 9, 186, 3, 3, 3,
    2, 5, 8, 3, 3, 3, 3, 3,
    0, 3, 0, 4, 6, 7, 9, 3,
    10, 3, 11, 2, 5, 8, 3, 3,
    3, 3, 3, 0, 2, 0, 2, 0,
    2, 0, 5, 0, 5, 0, 5, 0,
    4, 6, 7, 9, 3, 3, 3, 3,
    2, 5, 8, 3, 3, 3, 3, 3,
    0, 13, 15, 16, 18, 21, 12, 14,
    17, 196, 196, 196, 0, 196, 0, 12,
    0, 12, 0, 12, 0, 14, 0, 14,
    0, 14, 0, 13, 15, 16, 18, 19,
    12, 14, 17, 196, 196, 196, 0, 13,
    15, 16, 18, 12, 14, 17, 196, 196,
    196, 0, 22, 26, 44, 46, 48, 45,
    23, 23, 0, 22, 23, 23, 23, 0,
    22, 24, 23, 23, 23, 0, 25, 25,
    0, 197, 0, 22, 27, 24, 23, 23,
    23, 0, 28, 40, 42, 41, 0, 29,
    0, 30, 36, 38, 37, 0, 31, 0,
    25, 32, 34, 33, 0, 197, 33, 0,
    197, 25, 0, 35, 197, 33, 25, 0,
    197, 25, 0, 31, 37, 0, 31, 30,
    0, 31, 39, 37, 30, 0, 31, 30,
    0, 29, 41, 0, 29, 28, 0, 29,
    43, 41, 28, 0, 29, 28, 0, 22,
    27, 24, 45, 23, 23, 0, 22, 27,
    24, 26, 23, 23, 0, 22, 27, 47,
    24, 45, 26, 23, 23, 0, 22, 27,
    24, 26, 23, 23, 23, 0, 22, 24,
    49, 23, 23, 23, 0, 22, 24, 50,
    23, 23, 23, 0, 22, 51, 24, 23,
    23, 23, 0, 22, 52, 23, 23, 23,
    0, 185, 25, 53, 25, 53, 25, 25,
    53, 25, 0, 57, 197, 54, 54, 54,
    0, 57, 55, 55, 55, 0, 57, 56,
    56, 56, 0, 57, 0, 124, 58, 58,
    58, 0, 62, 59, 59, 59, 0, 62,
    60, 60, 60, 0, 62, 61, 61, 61,
    0, 62, 0, 124, 63, 63, 63, 0,
    67, 64, 64, 64, 0, 67, 65, 65,
    65, 0, 67, 66, 66, 66, 0, 67,
    0, 124, 68, 68, 68, 0, 72, 69,
    69, 69, 0, 72, 70, 70, 70, 0,
    72, 71, 71, 71, 0, 72, 0, 124,
    73, 73, 73, 0, 77, 74, 74, 74,
    0, 77, 75, 75, 75, 0, 77, 76,
    76, 76, 0, 77, 0, 98, 78, 78,
    78, 0, 82, 79, 79, 79, 0, 82,
    80, 80, 80, 0, 82, 81, 81, 81,
    0, 82, 0, 83, 91, 94, 98, 97,
    123, 123, 0, 27, 87, 84, 84, 84,
    0, 87, 85, 85, 85, 0, 87, 86,
    86, 86, 0, 87, 0, 88, 88, 88,
    0, 197, 89, 89, 89, 0, 197, 90,
    90, 90, 0, 197, 25, 25, 25, 0,
    27, 87, 92, 84, 84, 0, 27, 87,
    93, 85, 85, 0, 27, 87, 86, 86,
    86, 0, 27, 95, 87, 92, 96, 84,
    84, 0, 27, 87, 93, 85, 85, 85,
    0, 27, 87, 85, 85, 85, 0, 27,
    87, 96, 84, 84, 0, 197, 99, 99,
    99, 0, 103, 197, 100, 100, 100, 0,
    103, 197, 101, 101, 101, 0, 103, 197,
    102, 102, 102, 0, 103, 197, 0, 104,
    104, 104, 0, 108, 197, 105, 105, 105,
    0, 108, 197, 106, 106, 106, 0, 108,
    197, 107, 107, 107, 0, 108, 197, 0,
    109, 109, 109, 0, 113, 197, 110, 110,
    110, 0, 113, 197, 111, 111, 111, 0,
    113, 197, 112, 112, 112, 0, 113, 197,
    0, 114, 114, 114, 0, 118, 197, 115,
    115, 115, 0, 118, 197, 116, 116, 116,
    0, 118, 197, 117, 117, 117, 0, 118,
    197, 0, 119, 119, 119, 0, 87, 197,
    120, 120, 120, 0, 87, 197, 121, 121,
    121, 0, 87, 197, 122, 122, 122, 0,
    87, 197, 0, 87, 84, 84, 84, 0,
    125, 177, 180, 197, 183, 184, 184, 0,
    27, 129, 197, 126, 126, 126, 0, 129,
    197, 127, 127, 127, 0, 129, 197, 128,
    128, 128, 0, 129, 197, 0, 130, 169,
    172, 175, 176, 176, 0, 27, 134, 197,
    131, 131, 131, 0, 134, 197, 132, 132,
    132, 0, 134, 197, 133, 133, 133, 0,
    134, 197, 0, 135, 161, 164, 167, 168,
    168, 0, 27, 139, 197, 136, 136, 136,
    0, 139, 197, 137, 137, 137, 0, 139,
    197, 138, 138, 138, 0, 139, 197, 0,
    140, 153, 156, 159, 160, 160, 0, 27,
    144, 197, 141, 141, 141, 0, 144, 197,
    142, 142, 142, 0, 144, 197, 143, 143,
    143, 0, 144, 197, 0, 145, 146, 149,
    152, 119, 119, 0, 27, 87, 197, 120,
    120, 120, 0, 27, 87, 197, 147, 120,
    120, 0, 27, 87, 197, 148, 121, 121,
    0, 27, 87, 197, 122, 122, 122, 0,
    27, 150, 87, 197, 147, 151, 120, 120,
    0, 27, 87, 197, 148, 121, 121, 121,
    0, 27, 87, 197, 121, 121, 121, 0,
    27, 87, 197, 151, 120, 120, 0, 27,
    144, 197, 154, 141, 141, 0, 27, 144,
    197, 155, 142, 142, 0, 27, 144, 197,
    143, 143, 143, 0, 27, 157, 144, 197,
    154, 158, 141, 141, 0, 27, 144, 197,
    155, 142, 142, 142, 0, 27, 144, 197,
    142, 142, 142, 0, 27, 144, 197, 158,
    141, 141, 0, 144, 197, 141, 141, 141,
    0, 27, 139, 197, 162, 136, 136, 0,
    27, 139, 197, 163, 137, 137, 0, 27,
    139, 197, 138, 138, 138, 0, 27, 165,
    139, 197, 162, 166, 136, 136, 0, 27,
    139, 197, 163, 137, 137, 137, 0, 27,
    139, 197, 137, 137, 137, 0, 27, 139,
    197, 166, 136, 136, 0, 139, 197, 136,
    136, 136, 0, 27, 134, 197, 170, 131,
    131, 0, 27, 134, 197, 171, 132, 132,
    0, 27, 134, 197, 133, 133, 133, 0,
    27, 173, 134, 197, 170, 174, 131, 131,
    0, 27, 134, 197, 171, 132, 132, 132,
    0, 27, 134, 197, 132, 132, 132, 0,
    27, 134, 197, 174, 131, 131, 0, 134,
    197, 131, 131, 131, 0, 27, 129, 197,
    178, 126, 126, 0, 27, 129, 197, 179,
    127, 127, 0, 27, 129, 197, 128, 128,
    128, 0, 27, 181, 129, 197, 178, 182,
    126, 126, 0, 27, 129, 197, 179, 127,
    127, 127, 0, 27, 129, 197, 127, 127,
    127, 0, 27, 129, 197, 182, 126, 126,
    0, 129, 197, 126, 126, 126, 0, 124,
    197, 0, 188, 190, 191, 193, 194, 195,
    187, 189, 192, 186, 0, 186, 0, 187,
    0, 187, 0, 187, 0, 189, 0, 189,
    0, 189, 0, 11, 0, 186, 0, 13,
    15, 16, 18, 19, 20, 12, 14, 17,
    196, 196, 196, 0, 0, 0, 1, 2,
    3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26,
    27, 28, 29, 30, 31, 32, 33, 34,
    35, 36, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58,
    59, 60, 61, 62, 63, 64, 65, 66,
    67, 68, 69, 70, 71, 72, 73, 74,
    75, 76, 77, 78, 79, 80, 81, 82,
    83, 84, 85, 86, 87, 88, 89, 90,
    91, 92, 93, 94, 95, 96, 97, 98,
    99, 100, 101, 102, 103, 104, 105, 106,
    107, 108, 109, 110, 111, 112, 113, 114,
    115, 116, 117, 118, 119, 120, 121, 122,
    123, 124, 125, 126, 127, 128, 129, 130,
    131, 132, 133, 134, 135, 136, 137, 138,
    139, 140, 141, 142, 143, 144, 145, 146,
    147, 148, 149, 150, 151, 152, 153, 154,
    155, 156, 157, 158, 159, 160, 161, 162,
    163, 164, 165, 166, 167, 168, 169, 170,
    171, 172, 173, 174, 175, 176, 177, 178,
    179, 180, 181, 182, 183, 184, 185, 186,
    187, 188, 189, 190, 191, 192, 193, 194,
    195, 196, 197, 0};

static const signed char _address_cond_actions[] = {
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    3, 0, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 3, 0, 3, 0, 3,
    0, 3, 0, 3, 0, 3, 0, 3,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    3, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 3, 1, 3, 0,
    3, 0, 3, 0, 3, 0, 3, 0,
    3, 0, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 3, 0,
    0, 0, 0, 0, 0, 0, 1, 1,
    1, 3, 0, 0, 0, 0, 0, 0,
    0, 0, 3, 0, 0, 0, 0, 3,
    0, 0, 0, 0, 0, 3, 0, 0,
    3, 1, 3, 0, 0, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 3, 0,
    3, 0, 0, 0, 0, 3, 0, 3,
    0, 0, 0, 0, 3, 1, 0, 3,
    1, 0, 3, 0, 1, 0, 0, 3,
    1, 0, 3, 0, 0, 3, 0, 0,
    3, 0, 0, 0, 0, 3, 0, 0,
    3, 0, 0, 3, 0, 0, 3, 0,
    0, 0, 0, 3, 0, 0, 3, 0,
    0, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0,
    0, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    0, 0, 3, 0, 0, 0, 0, 0,
    3, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 3, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 3, 0, 3, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 3, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    3, 0, 3, 0, 0, 0, 0, 3,
    0, 0, 0, 0, 3, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 3, 0,
    3, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 3, 0, 0, 0, 0, 3,
    0, 0, 0, 0, 3, 0, 3, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    3, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 3, 0, 3, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 3, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    3, 0, 3, 0, 0, 0, 0, 0,
    0, 0, 3, 0, 0, 0, 0, 0,
    3, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 3, 0, 3, 0, 0, 0,
    3, 1, 0, 0, 0, 3, 1, 0,
    0, 0, 3, 1, 0, 0, 0, 3,
    0, 0, 0, 0, 0, 3, 0, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 0, 0, 0, 0, 0,
    3, 0, 0, 0, 0, 0, 3, 0,
    0, 0, 0, 0, 3, 1, 0, 0,
    0, 3, 0, 1, 0, 0, 0, 3,
    0, 1, 0, 0, 0, 3, 0, 1,
    0, 0, 0, 3, 0, 1, 3, 0,
    0, 0, 3, 0, 1, 0, 0, 0,
    3, 0, 1, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 3, 0, 1, 3,
    0, 0, 0, 3, 0, 1, 0, 0,
    0, 3, 0, 1, 0, 0, 0, 3,
    0, 1, 0, 0, 0, 3, 0, 1,
    3, 0, 0, 0, 3, 0, 1, 0,
    0, 0, 3, 0, 1, 0, 0, 0,
    3, 0, 1, 0, 0, 0, 3, 0,
    1, 3, 0, 0, 0, 3, 0, 1,
    0, 0, 0, 3, 0, 1, 0, 0,
    0, 3, 0, 1, 0, 0, 0, 3,
    0, 1, 3, 0, 0, 0, 0, 3,
    0, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 1, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 3, 0, 1, 0,
    0, 0, 3, 0, 1, 3, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 1,
    0, 0, 0, 3, 0, 1, 0, 0,
    0, 3, 0, 1, 0, 0, 0, 3,
    0, 1, 3, 0, 0, 0, 0, 0,
    0, 3, 0, 0, 1, 0, 0, 0,
    3, 0, 1, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 3, 0, 1, 3,
    0, 0, 0, 0, 0, 0, 3, 0,
    0, 1, 0, 0, 0, 3, 0, 1,
    0, 0, 0, 3, 0, 1, 0, 0,
    0, 3, 0, 1, 3, 0, 0, 0,
    0, 0, 0, 3, 0, 0, 1, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 3, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 0, 1, 0, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 1, 0, 0, 0, 3, 0,
    0, 1, 0, 0, 0, 3, 0, 0,
    1, 0, 0, 0, 3, 0, 0, 1,
    0, 0, 0, 3, 0, 0, 0, 1,
    0, 0, 0, 0, 3, 0, 0, 1,
    0, 0, 0, 0, 3, 0, 0, 1,
    0, 0, 0, 3, 0, 0, 1, 0,
    0, 0, 3, 0, 1, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 1, 0, 0, 0, 3, 0,
    0, 1, 0, 0, 0, 3, 0, 0,
    0, 1, 0, 0, 0, 0, 3, 0,
    0, 1, 0, 0, 0, 0, 3, 0,
    0, 1, 0, 0, 0, 3, 0, 0,
    1, 0, 0, 0, 3, 0, 1, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 3, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 0, 1, 0, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 0,
    3, 0, 0, 1, 0, 0, 0, 3,
    0, 0, 1, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 3, 0, 0, 1,
    0, 0, 0, 3, 0, 0, 1, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 3, 0, 0, 0, 1, 0, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 0, 3, 0, 0, 1, 0, 0,
    0, 3, 0, 0, 1, 0, 0, 0,
    3, 0, 1, 0, 0, 0, 3, 0,
    1, 3, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 3, 0,
    3, 0, 3, 0, 3, 0, 3, 0,
    3, 0, 3, 0, 3, 0, 3, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 3, 3, 0, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3,
    3, 0, 0, 0};

static const short _address_eof_trans[] = {
    1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093,
    1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101,
    1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109,
    1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117,
    1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125,
    1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133,
    1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141,
    1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149,
    1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157,
    1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165,
    1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173,
    1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181,
    1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189,
    1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197,
    1198, 1199, 1200, 1201, 1202, 1203, 1204, 1205,
    1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213,
    1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221,
    1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229,
    1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237,
    1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245,
    1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
    1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261,
    1262, 1263, 1264, 1265, 1266, 1267, 1268, 1269,
    1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277,
    1278, 1279, 1280, 1281, 1282, 1283, 0};

static const int address_start = 1;

bool is_address(const char *p, const char *pe)
{
	int cs = 0;

	const char *eof = pe;

	bool result = false;

	{
		cs = (int) address_start;
	}
	{
		int _klen;
		unsigned int _trans = 0;
		const signed char *_keys;
		const signed char *_acts;
		unsigned int _nacts;
	_resume : {
	}
		if (p == pe && p != eof)
			goto _out;
		if (p == eof) {
			if (_address_eof_trans[cs] > 0) {
				_trans = (unsigned int) _address_eof_trans[cs] - 1;
			}
		} else {
			_keys = (_address_trans_keys + (_address_key_offsets[cs]));
			_trans = (unsigned int) _address_index_offsets[cs];

			_klen = (int) _address_single_lengths[cs];
			if (_klen > 0) {
				const signed char *_lower = _keys;
				const signed char *_upper = _keys + _klen - 1;
				const signed char *_mid;
				while (1) {
					if (_upper < _lower) {
						_keys += _klen;
						_trans += (unsigned int) _klen;
						break;
					}

					_mid = _lower + ((_upper - _lower) >> 1);
					if (((*(p))) < (*(_mid)))
						_upper = _mid - 1;
					else if (((*(p))) > (*(_mid)))
						_lower = _mid + 1;
					else {
						_trans += (unsigned int) (_mid - _keys);
						goto _match;
					}
				}
			}

			_klen = (int) _address_range_lengths[cs];
			if (_klen > 0) {
				const signed char *_lower = _keys;
				const signed char *_upper = _keys + (_klen << 1) - 2;
				const signed char *_mid;
				while (1) {
					if (_upper < _lower) {
						_trans += (unsigned int) _klen;
						break;
					}

					_mid = _lower + (((_upper - _lower) >> 1) & ~1);
					if (((*(p))) < (*(_mid)))
						_upper = _mid - 2;
					else if (((*(p))) > (*(_mid + 1)))
						_lower = _mid + 2;
					else {
						_trans += (unsigned int) ((_mid - _keys) >> 1);
						break;
					}
				}
			}

		_match : {
		}
		}
		cs = (int) _address_cond_targs[_trans];

		if (_address_cond_actions[_trans] != 0) {

			_acts = (_address_actions + (_address_cond_actions[_trans]));
			_nacts = (unsigned int) (*(_acts));
			_acts += 1;
			while (_nacts > 0) {
				switch ((*(_acts))) {
				case 0: {
					{
						result = true;
					}
					break;
				}
				case 1: {
					{
						result = false;
					}
					break;
				}
				}
				_nacts -= 1;
				_acts += 1;
			}
		}

		if (p == eof) {
			if (cs >= 196)
				goto _out;
		} else {
			if (cs != 0) {
				p += 1;
				goto _resume;
			}
		}
	_out : {
	}
	}
	return result;
}
// #include <nlohmann/json-schema.hpp>


// #include "smtp-address-validator.hpp"


#include <algorithm>
#include <exception>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef JSON_SCHEMA_BOOST_REGEX
#	include <boost/regex.hpp>
#	define REGEX_NAMESPACE boost
#elif defined(JSON_SCHEMA_NO_REGEX)
#	define NO_STD_REGEX
#else
#	include <regex>
#	define REGEX_NAMESPACE std
#endif

/**
 * Many of the RegExes are from @see http://jmrware.com/articles/2009/uri_regexp/URI_regex.html
 */

namespace
{
template <typename T>
void range_check(const T value, const T min, const T max)
{
	if (!((value >= min) && (value <= max))) {
		std::stringstream out;
		out << "Value " << value << " should be in interval [" << min << "," << max << "] but is not!";
		throw std::invalid_argument(out.str());
	}
}

/** @see date_time_check */
void rfc3339_date_check(const std::string &value)
{
	const static REGEX_NAMESPACE::regex dateRegex{R"(^([0-9]{4})\-([0-9]{2})\-([0-9]{2})$)"};

	REGEX_NAMESPACE::smatch matches;
	if (!REGEX_NAMESPACE::regex_match(value, matches, dateRegex)) {
		throw std::invalid_argument(value + " is not a date string according to RFC 3339.");
	}

	const auto year = std::stoi(matches[1].str());
	const auto month = std::stoi(matches[2].str());
	const auto mday = std::stoi(matches[3].str());

	const auto isLeapYear = (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));

	range_check(month, 1, 12);
	if (month == 2) {
		range_check(mday, 1, isLeapYear ? 29 : 28);
	} else if (month <= 7) {
		range_check(mday, 1, month % 2 == 0 ? 30 : 31);
	} else {
		range_check(mday, 1, month % 2 == 0 ? 31 : 30);
	}
}

/** @see date_time_check */
void rfc3339_time_check(const std::string &value)
{
	const static REGEX_NAMESPACE::regex timeRegex{R"(^([0-9]{2})\:([0-9]{2})\:([0-9]{2})(\.[0-9]+)?(?:[Zz]|((?:\+|\-)[0-9]{2})\:([0-9]{2}))$)"};

	REGEX_NAMESPACE::smatch matches;
	if (!REGEX_NAMESPACE::regex_match(value, matches, timeRegex)) {
		throw std::invalid_argument(value + " is not a time string according to RFC 3339.");
	}

	auto hour = std::stoi(matches[1].str());
	auto minute = std::stoi(matches[2].str());
	auto second = std::stoi(matches[3].str());
	// const auto secfrac      = std::stof( matches[4].str() );

	range_check(hour, 0, 23);
	range_check(minute, 0, 59);

	int offsetHour = 0,
	    offsetMinute = 0;

	/* don't check the numerical offset if time zone is specified as 'Z' */
	if (!matches[5].str().empty()) {
		offsetHour = std::stoi(matches[5].str());
		offsetMinute = std::stoi(matches[6].str());

		range_check(offsetHour, -23, 23);
		range_check(offsetMinute, 0, 59);
		if (offsetHour < 0)
			offsetMinute *= -1;
	}

	/**
	 * @todo Could be made more exact by querying a leap second database and choosing the
	 *       correct maximum in {58,59,60}. This current solution might match some invalid dates
	 *       but it won't lead to false negatives. This only works if we know the full date, however
	 */

	auto day_minutes = hour * 60 + minute - (offsetHour * 60 + offsetMinute);
	if (day_minutes < 0)
		day_minutes += 60 * 24;
	hour = day_minutes % 24;
	minute = day_minutes / 24;

	if (hour == 23 && minute == 59)
		range_check(second, 0, 60); // possible leap-second
	else
		range_check(second, 0, 59);
}

/**
 * @see https://tools.ietf.org/html/rfc3339#section-5.6
 *
 * @verbatim
 * date-fullyear   = 4DIGIT
 * date-month      = 2DIGIT  ; 01-12
 * date-mday       = 2DIGIT  ; 01-28, 01-29, 01-30, 01-31 based on
 *                          ; month/year
 * time-hour       = 2DIGIT  ; 00-23
 * time-minute     = 2DIGIT  ; 00-59
 * time-second     = 2DIGIT  ; 00-58, 00-59, 00-60 based on leap second
 *                          ; rules
 * time-secfrac    = "." 1*DIGIT
 * time-numoffset  = ("+" / "-") time-hour ":" time-minute
 * time-offset     = "Z" / time-numoffset
 *
 * partial-time    = time-hour ":" time-minute ":" time-second
 *                  [time-secfrac]
 * full-date       = date-fullyear "-" date-month "-" date-mday
 * full-time       = partial-time time-offset
 *
 * date-time       = full-date "T" full-time
 * @endverbatim
 * NOTE: Per [ABNF] and ISO8601, the "T" and "Z" characters in this
 *       syntax may alternatively be lower case "t" or "z" respectively.
 */
void rfc3339_date_time_check(const std::string &value)
{
	const static REGEX_NAMESPACE::regex dateTimeRegex{R"(^([0-9]{4}\-[0-9]{2}\-[0-9]{2})[Tt]([0-9]{2}\:[0-9]{2}\:[0-9]{2}(?:\.[0-9]+)?(?:[Zz]|(?:\+|\-)[0-9]{2}\:[0-9]{2}))$)"};

	REGEX_NAMESPACE::smatch matches;
	if (!REGEX_NAMESPACE::regex_match(value, matches, dateTimeRegex)) {
		throw std::invalid_argument(value + " is not a date-time string according to RFC 3339.");
	}

	rfc3339_date_check(matches[1].str());
	rfc3339_time_check(matches[2].str());
}

const std::string decOctet{R"((?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9]))"}; // matches numbers 0-255
const std::string ipv4Address{"(?:" + decOctet + R"(\.){3})" + decOctet};
const std::string h16{R"([0-9A-Fa-f]{1,4})"};
const std::string h16Left{"(?:" + h16 + ":)"};
const std::string ipv6Address{
    "(?:"
    "(?:" +
    h16Left + "{6}"
              "|::" +
    h16Left + "{5}"
              "|(?:" +
    h16 + ")?::" + h16Left + "{4}"
                             "|(?:" +
    h16Left + "{0,1}" + h16 + ")?::" + h16Left + "{3}"
                                                 "|(?:" +
    h16Left + "{0,2}" + h16 + ")?::" + h16Left + "{2}"
                                                 "|(?:" +
    h16Left + "{0,3}" + h16 + ")?::" + h16Left +
    "|(?:" + h16Left + "{0,4}" + h16 + ")?::"
                                       ")(?:" +
    h16Left + h16 + "|" + ipv4Address + ")"
                                        "|(?:" +
    h16Left + "{0,5}" + h16 + ")?::" + h16 +
    "|(?:" + h16Left + "{0,6}" + h16 + ")?::"
                                       ")"};
const std::string ipvFuture{R"([Vv][0-9A-Fa-f]+\.[A-Za-z0-9\-._~!$&'()*+,;=:]+)"};
const std::string regName{R"((?:[A-Za-z0-9\-._~!$&'()*+,;=]|%[0-9A-Fa-f]{2})*)"};
const std::string host{
    "(?:"
    R"(\[(?:)" +
    ipv6Address + "|" + ipvFuture + R"()\])" +
    "|" + ipv4Address +
    "|" + regName +
    ")"};

const std::string uuid{R"([0-9a-fA-F]{8}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{12})"};

// from http://stackoverflow.com/questions/106179/regular-expression-to-match-dns-hostname-or-ip-address
const std::string hostname{R"(^([a-zA-Z0-9]|[a-zA-Z0-9][a-zA-Z0-9\-]{0,61}[a-zA-Z0-9])(\.([a-zA-Z0-9]|[a-zA-Z0-9][a-zA-Z0-9\-]{0,61}[a-zA-Z0-9]))*$)"};

bool is_ascii(std::string const &value)
{
	for (auto ch : value) {
		if (ch & 0x80) {
			return false;
		}
	}
	return true;
}

/**
 * @see
 *
 * @verbatim
 * URI           = scheme ":" hier-part [ "?" query ] [ "#" fragment ]
 *
 *  hier-part     = "//" authority path-abempty
 *               / path-absolute
 *               / path-rootless
 *               / path-empty
 *
 * URI-reference = URI / relative-ref
 *
 * absolute-URI  = scheme ":" hier-part [ "?" query ]
 *
 * relative-ref  = relative-part [ "?" query ] [ "#" fragment ]
 *
 * relative-part = "//" authority path-abempty
 *               / path-absolute
 *               / path-noscheme
 *               / path-empty
 *
 * scheme        = ALPHA *( ALPHA / DIGIT / "+" / "-" / "." )
 *
 * authority     = [ userinfo "@" ] host [ ":" port ]
 * userinfo      = *( unreserved / pct-encoded / sub-delims / ":" )
 * host          = IP-literal / IPv4address / reg-name
 * port          = *DIGIT
 *
 * IP-literal    = "[" ( IPv6address / IPvFuture  ) "]"
 *
 * IPvFuture     = "v" 1*HEXDIG "." 1*( unreserved / sub-delims / ":" )
 *
 * IPv6address   =                            6( h16 ":" ) ls32
 *               /                       "::" 5( h16 ":" ) ls32
 *               / [               h16 ] "::" 4( h16 ":" ) ls32
 *               / [ *1( h16 ":" ) h16 ] "::" 3( h16 ":" ) ls32
 *               / [ *2( h16 ":" ) h16 ] "::" 2( h16 ":" ) ls32
 *               / [ *3( h16 ":" ) h16 ] "::"    h16 ":"   ls32
 *               / [ *4( h16 ":" ) h16 ] "::"              ls32
 *               / [ *5( h16 ":" ) h16 ] "::"              h16
 *               / [ *6( h16 ":" ) h16 ] "::"
 *
 * h16           = 1*4HEXDIG
 * ls32          = ( h16 ":" h16 ) / IPv4address
 * IPv4address   = dec-octet "." dec-octet "." dec-octet "." dec-octet
 *    dec-octet     = DIGIT                 ; 0-9
 *               / %x31-39 DIGIT         ; 10-99
 *               / "1" 2DIGIT            ; 100-199
 *               / "2" %x30-34 DIGIT     ; 200-249
 *               / "25" %x30-35          ; 250-255
 *
 * reg-name      = *( unreserved / pct-encoded / sub-delims )
 *
 * path          = path-abempty    ; begins with "/" or is empty
 *               / path-absolute   ; begins with "/" but not "//"
 *               / path-noscheme   ; begins with a non-colon segment
 *               / path-rootless   ; begins with a segment
 *               / path-empty      ; zero characters
 *
 * path-abempty  = *( "/" segment )
 * path-absolute = "/" [ segment-nz *( "/" segment ) ]
 * path-noscheme = segment-nz-nc *( "/" segment )
 * path-rootless = segment-nz *( "/" segment )
 * path-empty    = 0<pchar>
 *
 * segment       = *pchar
 * segment-nz    = 1*pchar
 * segment-nz-nc = 1*( unreserved / pct-encoded / sub-delims / "@" )
 *               ; non-zero-length segment without any colon ":"
 *
 * pchar         = unreserved / pct-encoded / sub-delims / ":" / "@"
 *
 * query         = *( pchar / "/" / "?" )
 *
 * fragment      = *( pchar / "/" / "?" )
 *
 * pct-encoded   = "%" HEXDIG HEXDIG
 *
 * unreserved    = ALPHA / DIGIT / "-" / "." / "_" / "~"
 * reserved      = gen-delims / sub-delims
 * gen-delims    = ":" / "/" / "?" / "#" / "[" / "]" / "@"
 * sub-delims    = "!" / "$" / "&" / "'" / "(" / ")"
 *               / "*" / "+" / "," / ";" / "="
 *
 * @endverbatim
 * @see adapted from: https://github.com/jhermsmeier/uri.regex/blob/master/uri.regex
 *
 */
void rfc3986_uri_check(const std::string &value)
{
	const static std::string scheme{R"(([A-Za-z][A-Za-z0-9+\-.]*):)"};
	const static std::string hierPart{
	    R"((?:(\/\/)(?:((?:[A-Za-z0-9\-._~!$&'()*+,;=:]|)"
	    R"(%[0-9A-Fa-f]{2})*)@)?((?:\[(?:(?:(?:(?:[0-9A-Fa-f]{1,4}:){6}|)"
	    R"(::(?:[0-9A-Fa-f]{1,4}:){5}|)"
	    R"((?:[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){4}|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,1}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){3}|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,2}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){2}|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,3}[0-9A-Fa-f]{1,4})?::[0-9A-Fa-f]{1,4}:|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,4}[0-9A-Fa-f]{1,4})?::)(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|)"
	    R"((?:(?:25[0-5]|2[0-4][0-9]|)"
	    R"([01]?[0-9][0-9]?)\.){3}(?:25[0-5]|)"
	    R"(2[0-4][0-9]|)"
	    R"([01]?[0-9][0-9]?))|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,5}[0-9A-Fa-f]{1,4})?::[0-9A-Fa-f]{1,4}|)"
	    R"((?:(?:[0-9A-Fa-f]{1,4}:){0,6}[0-9A-Fa-f]{1,4})?::)|)"
	    R"([Vv][0-9A-Fa-f]+\.[A-Za-z0-9\-._~!$&'()*+,;=:]+)\]|)"
	    R"((?:(?:25[0-5]|)"
	    R"(2[0-4][0-9]|)"
	    R"([01]?[0-9][0-9]?)\.){3}(?:25[0-5]|)"
	    R"(2[0-4][0-9]|)"
	    R"([01]?[0-9][0-9]?)|)"
	    R"((?:[A-Za-z0-9\-._~!$&'()*+,;=]|)"
	    R"(%[0-9A-Fa-f]{2})*))(?::([0-9]*))?((?:\/(?:[A-Za-z0-9\-._~!$&'()*+,;=:@]|)"
	    R"(%[0-9A-Fa-f]{2})*)*)|)"
	    R"(\/((?:(?:[A-Za-z0-9\-._~!$&'()*+,;=:@]|)"
	    R"(%[0-9A-Fa-f]{2})+(?:\/(?:[A-Za-z0-9\-._~!$&'()*+,;=:@]|)"
	    R"(%[0-9A-Fa-f]{2})*)*)?)|)"
	    R"(((?:[A-Za-z0-9\-._~!$&'()*+,;=:@]|)"
	    R"(%[0-9A-Fa-f]{2})+(?:\/(?:[A-Za-z0-9\-._~!$&'()*+,;=:@]|)"
	    R"(%[0-9A-Fa-f]{2})*)*)|))"};

	const static std::string query{R"((?:\?((?:[A-Za-z0-9\-._~!$&'()*+,;=:@\/?]|%[0-9A-Fa-f]{2})*))?)"};
	const static std::string fragment{
	    R"((?:\#((?:[A-Za-z0-9\-._~!$&'()*+,;=:@\/?]|%[0-9A-Fa-f]{2})*))?)"};
	const static std::string uriFormat{scheme + hierPart + query + fragment};

	const static REGEX_NAMESPACE::regex uriRegex{uriFormat};

	if (!REGEX_NAMESPACE::regex_match(value, uriRegex)) {
		throw std::invalid_argument(value + " is not a URI string according to RFC 3986.");
	}
}

} // namespace

namespace nlohmann
{
namespace json_schema
{
/**
 * Checks validity for built-ins by converting the definitions given as ABNF in the linked RFC from
 * @see https://json-schema.org/understanding-json-schema/reference/string.html#built-in-formats
 * into regular expressions using @see https://www.msweet.org/abnf/ and some manual editing.
 *
 * @see https://json-schema.org/latest/json-schema-validation.html
 */
void default_string_format_check(const std::string &format, const std::string &value)
{
	if (format == "date-time") {
		rfc3339_date_time_check(value);
	} else if (format == "date") {
		rfc3339_date_check(value);
	} else if (format == "time") {
		rfc3339_time_check(value);
	} else if (format == "uri") {
		rfc3986_uri_check(value);
	} else if (format == "email") {
		if (!is_ascii(value)) {
			throw std::invalid_argument(value + " contains non-ASCII values, not RFC 5321 compliant.");
		}
		if (!is_address(&*value.begin(), &*value.end())) {
			throw std::invalid_argument(value + " is not a valid email according to RFC 5321.");
		}
	} else if (format == "idn-email") {
		if (!is_address(&*value.begin(), &*value.end())) {
			throw std::invalid_argument(value + " is not a valid idn-email according to RFC 6531.");
		}
	} else if (format == "hostname") {
		static const REGEX_NAMESPACE::regex hostRegex{hostname};
		if (!REGEX_NAMESPACE::regex_match(value, hostRegex)) {
			throw std::invalid_argument(value + " is not a valid hostname according to RFC 3986 Appendix A.");
		}
	} else if (format == "ipv4") {
		const static REGEX_NAMESPACE::regex ipv4Regex{"^" + ipv4Address + "$"};
		if (!REGEX_NAMESPACE::regex_match(value, ipv4Regex)) {
			throw std::invalid_argument(value + " is not an IPv4 string according to RFC 2673.");
		}
	} else if (format == "ipv6") {
		static const REGEX_NAMESPACE::regex ipv6Regex{ipv6Address};
		if (!REGEX_NAMESPACE::regex_match(value, ipv6Regex)) {
			throw std::invalid_argument(value + " is not an IPv6 string according to RFC 5954.");
		}
	} else if (format == "uuid") {
		static const REGEX_NAMESPACE::regex uuidRegex{uuid};
		if (!REGEX_NAMESPACE::regex_match(value, uuidRegex)) {
			throw std::invalid_argument(value + " is not an uuid string according to RFC 4122.");
		}
	} else if (format == "regex") {
		try {
			REGEX_NAMESPACE::regex re(value, std::regex::ECMAScript);
		} catch (std::exception &exception) {
			throw exception;
		}
	} else {
		/* yet unsupported JSON schema draft 7 built-ins */
		static const std::vector<std::string> jsonSchemaStringFormatBuiltIns{
		    "date-time", "time", "date", "email", "idn-email", "hostname", "idn-hostname", "ipv4", "ipv6", "uri",
		    "uri-reference", "iri", "iri-reference", "uri-template", "json-pointer", "relative-json-pointer", "regex"};
		if (std::find(jsonSchemaStringFormatBuiltIns.begin(), jsonSchemaStringFormatBuiltIns.end(), format) != jsonSchemaStringFormatBuiltIns.end()) {
			throw std::logic_error("JSON schema string format built-in " + format + " not yet supported. " +
			                       "Please open an issue or use a custom format checker.");
		}

		throw std::logic_error("Don't know how to validate " + format);
	}
}
} // namespace json_schema
} // namespace nlohmann
