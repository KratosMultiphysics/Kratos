# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as HDF5

# --- STD Imports ---
import pathlib
import typing


def NoOpExtractor(model: KratosMultiphysics.Model) -> KratosMultiphysics.Parameters:
    return KratosMultiphysics.Parameters()


class RecordedPredicate:
    """ @brief Helper class for storing storing and recalling the results of a predicate in order."""

    def __init__(self, predicate: typing.Callable):
        self.__predicate = predicate
        self.__indices = []
        self.__index = 0
        self.__max_recorded_index = 0
        self.__recording = False

    def StartRecording(self) -> None:
        self.__recording = True
        self.__index = 0
        self.__indices = []
        self.__max_recorded_index = 0

    def StopRecording(self) -> None:
        self.__recording = False
        self.__index = 0

    def __call__(self, model: KratosMultiphysics.Model) -> bool:
        if self.__recording:
            result = self.__predicate(model)
            if result:
                self.__indices.append(self.__index)
            self.__max_recorded_index += 1
        else:
            if self.__max_recorded_index < self.__index:
                raise RuntimeError(f"Predicate record out of range {self.__index} ({self.__max_recorded_index})")
            result = self.__index in self.__indices

        self.__index += 1
        return result


class OutputJournal:
    """ @brief A @ref Journal for tracking and managing output files.
        @details @a OutputJournal shares most functionality with @ref Journal,
                 but each entry represents a file on the file system. Entries
                 must have a 'file_path' string item that points to its associated
                 file. Apart from having no line breaks, there are no further
                 restrictions on the res of the JSON object.
        @details Adding a entry does not create a new file, but erasing entries
                 or clearing the @a OutputJournal deletes the associated files
                 from the system.
        @warning If a new entry's 'file_path' matches that of an existing entry,
                 the existing entry is overwritten and a warning is issued.
        @throws If there is no 'file_path' in the resulting @ref Parameters object.
    """

    def __init__(self, journal_path: pathlib.Path, extractor: typing.Callable = NoOpExtractor):
        self.__journal = HDF5.Journal(journal_path, extractor)

        # Validate existing journal file
        entries = {}
        for entry in self.__journal:
            self.__ValidateEntryFormat(entry)
            file_path = entry["file_path"].GetString()
            entries.setdefault(file_path, [])
            entries[file_path].append(entry)

        for entry_list in entries.values():
            if 1 < len(entry_list):
                self.__PruneDuplicatesOf(entry[-1])


    def GetFilePath(self) -> pathlib.Path:
        return pathlib.Path(self.__journal.GetFilePath())


    def SetExtractor(self, extractor: typing.Callable) -> None:
        self.__journal.SetExtractor(extractor)


    def Push(self, model: KratosMultiphysics.Model) -> None:
        """ @brief Parse the input model and write the result to the journal, then check its format and uniqueness."""
        # Push a new entry without checking anything new.
        self.__journal.Push(model)

        # Check whether the entry that just got pushed is
        # of the expected format.
        last_entry = self.__GetLastEntry()
        self.__ValidateEntryFormat(last_entry)

        # Now check whether there are any duplicates,
        # and keep only the most recent ones if there are any.
        entries = {}
        require_erase = False

        for entry in self.__journal:
            path = entry["file_path"].GetString()
            if path in entries:
                require_erase = True
            entries[path] = entry

        if require_erase:
            self.__PruneDuplicatesOf(last_entry)


    def EraseIf(self, predicate: typing.Callable) -> None:
        """ @brief Remove all entries from the journal that match the input predicate."""
        # The input predicate may store state, so it must only be run once
        # ==> a new index-based predicate is constructed and passed instead
        # (this assumes that iteration over the journal happens in the same order)
        recorded_predicate = RecordedPredicate(predicate)
        recorded_predicate.StartRecording()

        for entry in self.__journal:
            if recorded_predicate(entry):
                path = pathlib.Path(entry["file_path"].GetString())
                if path.is_file():
                    path.unlink()
                elif path.is_dir():
                    raise FileExistsError(f"{path} is a directory")

        recorded_predicate.StopRecording()
        self.__journal.EraseIf(recorded_predicate)


    def Clear(self) -> None:
        """ @brief Remove the journal file and every file it points to.
            @throws If one of its entries points to a directory.
        """
        for entry in self.__journal:
            path = pathlib.Path(entry["file_path"].GetString())
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                raise FileExistsError(f"{path} is a directory")
        self.__journal.Clear()


    def __GetLastEntry(self) -> KratosMultiphysics.Parameters:
        # The reason why this check is so ugly is because Journal's
        # iterator is not random access.
        size = len(self)
        if size:
            last_index = size - 1
            for index, entry in enumerate(self):
                if index == last_index:
                    return entry


    def __PruneDuplicatesOf(self, entry_to_keep: KratosMultiphysics.Parameters) -> None:
        """ @brief Delete entries that match @a file_path except one that is identical to @a entry_to_keep."""
        class ErasePredicate:
            def __init__(self, entry_to_keep: KratosMultiphysics.Parameters):
                self.__found_exact_match = False
                self.__entry_to_keep = entry_to_keep
                self.__target_file_path = entry_to_keep["file_path"].GetString()

            def __call__(self, entry: KratosMultiphysics.Parameters) -> bool:
                if entry["file_path"].GetString() == self.__target_file_path:
                    if entry.IsEquivalentTo(self.__entry_to_keep):
                        if self.__found_exact_match:
                            return True
                        else:
                            self.__found_exact_match = True
                            return False
                    else:
                        return True
                else:
                    return False

        self.EraseIf(ErasePredicate(entry_to_keep))


    @staticmethod
    def __ValidateEntryFormat(entry: KratosMultiphysics.Parameters) -> None:
        if not entry.Has("file_path") or not entry["file_path"].IsString():
            raise ValueError(f"Invalid entry (missing or invalid 'file_path'): {entry}")


    def __len__(self) -> int:
        return self.__journal.__len__()


    def __iter__(self) -> typing.Iterable:
        return self.__journal.__iter__()
