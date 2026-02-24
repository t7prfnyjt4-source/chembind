import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      "/api": {
        target: "http://127.0.0.1:8000",
        changeOrigin: true,
        // IMPORTANT: do NOT rewrite /api away.
        // Backend expects /api/* routes.
      },
    },
  },
});
