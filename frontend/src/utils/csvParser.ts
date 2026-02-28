export type CsvRow = { rowId: string; smiles: string };

function splitCsvLine(line: string): string[] {
  const out: string[] = [];
  let cur = "";
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const ch = line[i];

    if (inQuotes) {
      if (ch === '"') {
        const next = line[i + 1];
        if (next === '"') {
          cur += '"';
          i++;
        } else {
          inQuotes = false;
        }
      } else {
        cur += ch;
      }
    } else {
      if (ch === ",") {
        out.push(cur);
        cur = "";
      } else if (ch === '"') {
        inQuotes = true;
      } else {
        cur += ch;
      }
    }
  }

  out.push(cur);
  return out.map((s) => s.trim());
}

function normalizeLines(csvText: string): string[] {
  const text = csvText
    .replace(/^\uFEFF/, "")
    .replace(/\r\n/g, "\n")
    .replace(/\r/g, "\n");
  return text.split("\n");
}

export function parseSmilesCsv(csvText: string): CsvRow[] {
  const lines = normalizeLines(csvText).map((l) => l.trim());
  const nonEmptyIdx = lines.findIndex((l) => l.length > 0);
  if (nonEmptyIdx === -1) return [];

  const firstLine = lines[nonEmptyIdx];
  const firstCells = splitCsvLine(firstLine).map((c) =>
    c.toLowerCase().trim()
  );

  let startRow = nonEmptyIdx;
  let smilesCol = 0;

  if (firstCells.includes("smiles")) {
    smilesCol = firstCells.indexOf("smiles");
    startRow = nonEmptyIdx + 1;
  }

  const rows: CsvRow[] = [];

  for (let i = startRow; i < lines.length; i++) {
    const raw = lines[i];
    if (!raw) continue;

    const cells = splitCsvLine(raw);
    const smiles = (cells[smilesCol] ?? "").trim();
    if (!smiles) continue;

    rows.push({
      rowId: `row_${i + 1}`,
      smiles,
    });
  }

  return rows;
}
