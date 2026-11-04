# hapapp-utils init

Write a commented settings file to fill in.

## Synopsis

```bash
hapapp-utils init [path] [--workflow {assign,build-db}] [--force]
```

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `path` | `run.yaml` | Where to write the file |
| `--workflow` | `assign` | Which command the settings are for: `assign` or `build-db` |
| `--force` | off | Replace the file if it already exists |
| `--debug` | off | Show the full Python error if something unexpected fails |

## Examples

```bash
hapapp-utils init run.yaml
hapapp-utils init build.yaml --workflow build-db
```

The file comments out every optional setting and shows its default. Remove the `#` from any setting you
want to change. [Settings](../reference/settings.md) describes each one.
