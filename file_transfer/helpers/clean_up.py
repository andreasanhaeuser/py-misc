def delete_empty_dirs(dirname):
    """Recursively delete empty directories."""
    dirname = dirname.rstrip('/')
    assert os.path.isdir(dirname)

    # check if directory is empty
    contents = os.listdir(dirname)
    if len(contents) > 0:
        return

    parent = os.path.dirname(dirname)

    # remove directory
    print('Remove empty directory: %s' % dirname)
    os.rmdir(dirname)

    # proceed with parent
    return delete_empty_dirs(parent)
