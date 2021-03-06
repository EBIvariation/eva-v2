/*
 * Copyright 2016 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *          http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package embl.ebi.variation.eva.pipeline.listener;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.batch.core.SkipListener;

/**
 * @author Diego Poggioli
 *
 * Just logs skipped items into a differnt file (different logger)
 */
public class SkipCheckingListener implements SkipListener {
    private static final Logger logger = LoggerFactory.getLogger("SkipLogging");

    @Override
    public void onSkipInRead(Throwable t) {
        logger.error("Skipped line during READ step: " + t.getMessage());
    }

    @Override
    public void onSkipInWrite(Object item, Throwable t) {
        logger.error("Skipped line during WRITE: " + t.getMessage());
    }

    @Override
    public void onSkipInProcess(Object item, Throwable t) {
        logger.error("Skipped line during PROCESS: " + t.getMessage());
    }
}
